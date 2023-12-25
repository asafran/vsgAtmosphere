#version 450
#extension GL_ARB_separate_shader_objects : enable
#pragma import_defines (ATMOSHPERE_RADIANCE, ATMOSHPERE_VIEWER_IN_SPACE, ATMOSHPERE_COMBINED_SCATTERING_TEXTURES, VSG_POINT_SPRITE, VSG_DIFFUSE_MAP, VSG_GREYSCALE_DIFFUSE_MAP, VSG_EMISSIVE_MAP, VSG_LIGHTMAP_MAP, VSG_NORMAL_MAP, VSG_SPECULAR_MAP, VSG_TWO_SIDED_LIGHTING, SHADOWMAP_DEBUG)

#include "rendering_constants.glsl"
#include "constants.glsl"
#include "functions.glsl"

#ifdef VSG_DIFFUSE_MAP
layout(set = MATERIAL_DESCRIPTOR_SET, binding = 0) uniform sampler2D diffuseMap;
#endif

#ifdef VSG_NORMAL_MAP
layout(set = MATERIAL_DESCRIPTOR_SET, binding = 2) uniform sampler2D normalMap;
#endif

#ifdef VSG_LIGHTMAP_MAP
layout(set = MATERIAL_DESCRIPTOR_SET, binding = 3) uniform sampler2D aoMap;
#endif

#ifdef VSG_EMISSIVE_MAP
layout(set = MATERIAL_DESCRIPTOR_SET, binding = 4) uniform sampler2D emissiveMap;
#endif

#ifdef VSG_SPECULAR_MAP
layout(set = MATERIAL_DESCRIPTOR_SET, binding = 5) uniform sampler2D specularMap;
#endif

layout(set = MATERIAL_DESCRIPTOR_SET, binding = 10) uniform MaterialData
{
    vec4 ambientColor;
    vec4 diffuseColor;
    vec4 specularColor;
    vec4 emissiveColor;
    float shininess;
    float alphaMask;
    float alphaMaskCutoff;
} material;

layout(set = VIEW_DESCRIPTOR_SET, binding = 0) uniform LightData
{
    vec4 values[2048];
} lightData;


layout(set = VIEW_DESCRIPTOR_SET, binding = 2) uniform sampler2DArrayShadow shadowMaps;

layout(location = 0) in vec3 eyePos;
layout(location = 1) in vec3 normalDir;
layout(location = 2) in vec4 vertexColor;
#ifndef VSG_POINT_SPRITE
layout(location = 3) in vec2 texCoord0;
#endif
layout(location = 5) in vec3 viewRay;

layout(location = 0) out vec4 outColor;

// Find the normal for this fragment, pulling either from a predefined normal map
// or from the interpolated mesh normal and tangent attributes.
vec3 getNormal()
{
    vec3 result;
#ifdef VSG_NORMAL_MAP
    // Perturb normal, see http://www.thetenthplanet.de/archives/1180
    vec3 tangentNormal = texture(normalMap, texCoord0).xyz * 2.0 - 1.0;

    //tangentNormal *= vec3(2,2,1);

    vec3 q1 = dFdx(eyePos);
    vec3 q2 = dFdy(eyePos);
    vec2 st1 = dFdx(texCoord0);
    vec2 st2 = dFdy(texCoord0);

    vec3 N = normalize(normalDir);
    vec3 T = normalize(q1 * st2.t - q2 * st1.t);
    vec3 B = -normalize(cross(N, T));
    mat3 TBN = mat3(T, B, N);

    result = normalize(TBN * tangentNormal);
#else
    result = normalize(normalDir);
#endif
#ifdef VSG_TWO_SIDED_LIGHTING
    if (!gl_FrontFacing)
        result = -result;
#endif
    return result;
}

layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 0) uniform sampler2D transmittance_texture;
layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 1) uniform sampler2D irradiance_texture;
layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 2) uniform sampler3D scattering_texture;
layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 3) uniform sampler3D single_mie_scattering_texture;

layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 4, std140) uniform Settings
{
	vec4 whitePoint;
	vec2 sunSize;
} settings;

layout(set = POSITIONAL_DESCRIPTOR_SET, binding = 0, std140) uniform Positional
{
	vec3 sunDirection;
    float exposure;
    vec3 cameraPos;
    float radius;
    float mu_s;
} positional;

#ifdef ATMOSHPERE_VIEWER_IN_SPACE
RadianceSpectrum GetSkyRadianceToPoint(IN(Position) view_ray, IN(Direction) view_direction, OUT(DimensionlessSpectrum) transmittance)
{
    // Compute the distance to the top atmosphere boundary along the view ray,
    // assuming the viewer is in space (or NaN if the view ray does not intersect
    // the atmosphere).
    Position camera = positional.cameraPos;
    Length r = positional.radius;
    Length rmu = dot(camera, view_direction);
    Length distance_to_top_atmosphere_boundary = -rmu -
        sqrt(rmu * rmu - r * r + top_radius * top_radius);
    // If the viewer is in space and the view ray intersects the  move
    // the viewer to the top atmosphere boundary (along the view ray):
    if (distance_to_top_atmosphere_boundary > 0.0 * m) {
        camera = camera + view_direction * distance_to_top_atmosphere_boundary;
        r = top_radius;
        rmu += distance_to_top_atmosphere_boundary;
    }

    // Compute the r, mu, mu_s and nu parameters for the first texture lookup.
    Number mu = rmu / r;
    Number mu_s = dot(camera, sun_direction) / r;
    Number nu = dot(view_direction, sun_direction);
    Length d = length(view_ray);
    bool ray_r_mu_intersects_ground = RayIntersectsGround( r, mu);

    transmittance = GetTransmittance( transmittance_texture,
        r, mu, d, ray_r_mu_intersects_ground);

    IrradianceSpectrum single_mie_scattering;
    IrradianceSpectrum scattering = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);

    // Compute the r, mu, mu_s and nu parameters for the second texture lookup.
    // If shadow_length is not 0 (case of light shafts), we want to ignore the
    // scattering along the last shadow_length meters of the view ray, which we
    // do by subtracting shadow_length from d (this way scattering_p is equal to
    // the S|x_s=x_0-lv term in Eq. (17) of our paper).
    Length r_p = ClampRadius( sqrt(d * d + 2.0 * r * mu * d + r * r));
    Number mu_p = (r * mu + d) / r_p;
    Number mu_s_p = (r * mu_s + d * nu) / r_p;

    IrradianceSpectrum single_mie_scattering_p;
    IrradianceSpectrum scattering_p = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
        single_mie_scattering_p);

    // Combine the lookup results to get the scattering between camera and point.
    scattering = scattering - transmittance * scattering_p;
    single_mie_scattering =
        single_mie_scattering - transmittance * single_mie_scattering_p;
    #ifdef ATMOSHPERE_COMBINED_SCATTERING_TEXTURES
    single_mie_scattering = GetExtrapolatedSingleMieScattering(
        vec4(scattering, single_mie_scattering.r));
    #endif

    // Hack to avoid rendering artifacts when the sun is below the horizon.
    single_mie_scattering = single_mie_scattering *
        smoothstep(Number(0.0), Number(0.01), mu_s);

    return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
        MiePhaseFunction(mie_phase_function_g, nu);
}
#else
RadianceSpectrum GetSkyRadianceToPoint( IN(Position) view_ray, IN(Direction) view_direction, OUT(DimensionlessSpectrum) transmittance)
{
    // Compute the distance to the top atmosphere boundary along the view ray,
    // assuming the viewer is in space (or NaN if the view ray does not intersect
    // the atmosphere).
    Length r = positional.radius;
    Length rmu = dot(positional.cameraPos, view_direction);

    // Compute the r, mu, mu_s and nu parameters for the first texture lookup.
    Number mu = rmu / r;
    //Number mu_s = positional.mu_s;
    Number mu_s = dot(positional.cameraPos, positional.sunDirection) / r;
    Number nu = dot(view_direction, positional.sunDirection);
    Length d = length(view_ray);
    bool ray_r_mu_intersects_ground = RayIntersectsGround( r, mu);

    transmittance = GetTransmittance( transmittance_texture,
        r, mu, d, ray_r_mu_intersects_ground);

    IrradianceSpectrum single_mie_scattering;
    IrradianceSpectrum scattering = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);

    // Compute the r, mu, mu_s and nu parameters for the second texture lookup.
    // If shadow_length is not 0 (case of light shafts), we want to ignore the
    // scattering along the last shadow_length meters of the view ray, which we
    // do by subtracting shadow_length from d (this way scattering_p is equal to
    // the S|x_s=x_0-lv term in Eq. (17) of our paper).
    Length r_p = ClampRadius( sqrt(d * d + 2.0 * r * mu * d + r * r));
    Number mu_p = (r * mu + d) / r_p;
    Number mu_s_p = (r * mu_s + d * nu) / r_p;

    IrradianceSpectrum single_mie_scattering_p;
    IrradianceSpectrum scattering_p = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
        single_mie_scattering_p);

    // Combine the lookup results to get the scattering between camera and point.
    scattering = scattering - transmittance * scattering_p;
    single_mie_scattering =
        single_mie_scattering - transmittance * single_mie_scattering_p;
    #ifdef ATMOSHPERE_COMBINED_SCATTERING_TEXTURES
    single_mie_scattering = GetExtrapolatedSingleMieScattering(
        vec4(scattering, single_mie_scattering.r));
    #endif

    // Hack to avoid rendering artifacts when the sun is below the horizon.
    single_mie_scattering = single_mie_scattering *
        smoothstep(Number(0.0), Number(0.01), mu_s);

    return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
        MiePhaseFunction(mie_phase_function_g, nu);
    }
    #endif

    IrradianceSpectrum GetSunAndSkyIrradiance( IN(Position) point, IN(Direction) normal, OUT(IrradianceSpectrum) sky_irradiance)
    {
    Length r = length(point);
    Number mu_s = dot(point, positional.sunDirection) / r;

    // Indirect irradiance (approximated if the surface is not horizontal).
    sky_irradiance = GetIrradiance( irradiance_texture, r, mu_s) *
        (1.0 + dot(normal, point) / r) * 0.5;

    // Direct irradiance.
    return solar_irradiance.xyz *
        GetTransmittanceToSun(
            transmittance_texture, r, mu_s) *
        max(dot(normal, positional.sunDirection), 0.0);
}

vec3 computeLighting(vec3 ambientColor, vec3 diffuseColor, vec3 specularColor, vec3 emissiveColor, float shininess, float ambientOcclusion, vec3 ld, vec3 nd, vec3 vd)
{
    vec3 color = vec3(0.0);
    color.rgb += ambientColor;

    float diff = max(dot(ld, nd), 0.0);
    color.rgb += diffuseColor * diff;

    if (diff > 0.0)
    {
        vec3 halfDir = normalize(ld + vd);
        color.rgb += specularColor * pow(max(dot(halfDir, nd), 0.0), shininess);
    }

    vec3 result = color + emissiveColor;
    result *= ambientOcclusion;

    return result;
}

void main()
{
    float brightnessCutoff = 0.001;

#ifdef VSG_POINT_SPRITE
    vec2 texCoord0 = gl_PointCoord.xy;
#endif

    vec4 diffuseColor = vertexColor * material.diffuseColor;
#ifdef VSG_DIFFUSE_MAP
    #ifdef VSG_GREYSCALE_DIFFUSE_MAP
        float v = texture(diffuseMap, texCoord0.st).s;
        diffuseColor *= vec4(v, v, v, 1.0);
    #else
        diffuseColor *= texture(diffuseMap, texCoord0.st);
    #endif
#endif

    vec4 ambientColor = diffuseColor * material.ambientColor * material.ambientColor.a;
    vec4 specularColor = material.specularColor;
    vec4 emissiveColor = material.emissiveColor;
    float shininess = material.shininess;
    float ambientOcclusion = 1.0;

    if (material.alphaMask == 1.0f)
    {
        if (diffuseColor.a < material.alphaMaskCutoff)
            discard;
    }

#ifdef VSG_EMISSIVE_MAP
    emissiveColor *= texture(emissiveMap, texCoord0.st);
#endif

#ifdef VSG_LIGHTMAP_MAP
    ambientOcclusion *= texture(aoMap, texCoord0.st).r;
#endif

#ifdef VSG_SPECULAR_MAP
    specularColor *= texture(specularMap, texCoord0.st);
#endif

    vec3 nd = getNormal();
    vec3 vd = normalize(viewRay);

    vec3 color = vec3(0.0, 0.0, 0.0);

    vec4 lightNums = lightData.values[0];
    int numAmbientLights = int(lightNums[0]);
    int numDirectionalLights = int(lightNums[1]);
    int numPointLights = int(lightNums[2]);
    int numSpotLights = int(lightNums[3]);
    int index = 1;

    if (numAmbientLights>0)
    {
        // ambient lights
        for(int i = 0; i<numAmbientLights; ++i)
        {
            vec4 lightColor = lightData.values[index++];
            color += (ambientColor.rgb * lightColor.rgb) * (lightColor.a);
        }
    }

    // index used to step through the shadowMaps array
    int shadowMapIndex = 0;

    if (numDirectionalLights>0)
    {
        // directional lights
        for(int i = 0; i<numDirectionalLights; ++i)
        {
            vec4 lightColor = lightData.values[index++];
            vec3 direction = -lightData.values[index++].xyz;
            vec4 shadowMapSettings = lightData.values[index++];

            float brightness = lightColor.a;

            // check shadow maps if required
            bool matched = false;
            while ((shadowMapSettings.r > 0.0 && brightness > brightnessCutoff) && !matched)
            {
                mat4 sm_matrix = mat4(lightData.values[index++],
                                      lightData.values[index++],
                                      lightData.values[index++],
                                      lightData.values[index++]);

                vec4 sm_tc = (sm_matrix) * vec4(eyePos, 1.0);

                if (sm_tc.x >= 0.0 && sm_tc.x <= 1.0 && sm_tc.y >= 0.0 && sm_tc.y <= 1.0 && sm_tc.z >= 0.0 /* && sm_tc.z <= 1.0*/)
                {
                    matched = true;

                    float coverage = texture(shadowMaps, vec4(sm_tc.st, shadowMapIndex, sm_tc.z)).r;
                    brightness *= (1.0-coverage);

#ifdef SHADOWMAP_DEBUG
                    if (shadowMapIndex==0) color = vec3(1.0, 0.0, 0.0);
                    else if (shadowMapIndex==1) color = vec3(0.0, 1.0, 0.0);
                    else if (shadowMapIndex==2) color = vec3(0.0, 0.0, 1.0);
                    else if (shadowMapIndex==3) color = vec3(1.0, 1.0, 0.0);
                    else if (shadowMapIndex==4) color = vec3(0.0, 1.0, 1.0);
                    else color = vec3(1.0, 1.0, 1.0);
#endif
                }

                ++shadowMapIndex;
                shadowMapSettings.r -= 1.0;
            }

            if (shadowMapSettings.r > 0.0)
            {
                // skip lightData and shadowMap entries for shadow maps that we haven't visited for this light
                // so subsequent light pointions are correct.
                index += 4 * int(shadowMapSettings.r);
                shadowMapIndex += int(shadowMapSettings.r);
            }

            // if light is too dim/shadowed to effect the rendering skip it
            if (brightness <= brightnessCutoff ) continue;

            float unclamped_LdotN = dot(direction, nd);

            float diff = max(unclamped_LdotN, 0.0);
            color.rgb += (diffuseColor.rgb * lightColor.rgb) * (diff * brightness);

            if (shininess > 0.0 && diff > 0.0)
            {
                vec3 halfDir = normalize(direction + vd);
                color.rgb += specularColor.rgb * (pow(max(dot(halfDir, nd), 0.0), shininess) * brightness);
            }
        }
    }

    if (numPointLights>0)
    {
        // point light
        for(int i = 0; i<numPointLights; ++i)
        {
            vec4 lightColor = lightData.values[index++];
            vec3 position = lightData.values[index++].xyz;

            float brightness = lightColor.a;

            // if light is too dim/shadowed to effect the rendering skip it
            if (brightness <= brightnessCutoff ) continue;

            vec3 delta = position - eyePos;
            float distance2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
            vec3 direction = delta / sqrt(distance2);
            float scale = brightness / distance2;

            float unclamped_LdotN = dot(direction, nd);

            float diff = scale * max(unclamped_LdotN, 0.0);

            color.rgb += (diffuseColor.rgb * lightColor.rgb) * diff;
            if (shininess > 0.0 && diff > 0.0)
            {
                vec3 halfDir = normalize(direction + vd);
                color.rgb += specularColor.rgb * (pow(max(dot(halfDir, nd), 0.0), shininess) * scale);
            }
        }
    }

    if (numSpotLights>0)
    {
        // spot light
        for(int i = 0; i<numSpotLights; ++i)
        {
            vec4 lightColor = lightData.values[index++];
            vec4 position_cosInnerAngle = lightData.values[index++];
            vec4 lightDirection_cosOuterAngle = lightData.values[index++];

            float brightness = lightColor.a;

            // if light is too dim/shadowed to effect the rendering skip it
            if (brightness <= brightnessCutoff ) continue;

            vec3 delta = position_cosInnerAngle.xyz - eyePos;
            float distance2 = delta.x * delta.x + delta.y * delta.y + delta.z * delta.z;
            vec3 direction = delta / sqrt(distance2);

            float dot_lightdirection = dot(lightDirection_cosOuterAngle.xyz, -direction);
            float scale = (brightness  * smoothstep(lightDirection_cosOuterAngle.w, position_cosInnerAngle.w, dot_lightdirection)) / distance2;

            float unclamped_LdotN = dot(direction, nd);

            float diff = scale * max(unclamped_LdotN, 0.0);
            color.rgb += (diffuseColor.rgb * lightColor.rgb) * diff;
            if (shininess > 0.0 && diff > 0.0)
            {
                vec3 halfDir = normalize(direction + vd);
                color.rgb += specularColor.rgb * (pow(max(dot(halfDir, nd), 0.0), shininess) * scale);
            }
        }
    }

    vec3 scaledViewRay = viewRay / 1000.0;

    vec3 point = positional.cameraPos + scaledViewRay;

    vec3 sky_irradiance;
	vec3 sun_irradiance = GetSunAndSkyIrradiance(point, -nd, sky_irradiance);
#ifdef ATMOSHPERE_RADIANCE
    sky_irradiance *= SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
    sun_irradiance *= SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
#endif
	vec3 sphere_radiance = (1.0 / PI) * (sun_irradiance + sky_irradiance);

	vec3 transmittance;
	vec3 in_scatter = GetSkyRadianceToPoint(scaledViewRay , vd, transmittance);
#ifdef ATMOSHPERE_RADIANCE
    in_scatter *= SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
#endif
	sphere_radiance = sphere_radiance * transmittance + in_scatter;

	sphere_radiance = pow(vec3(1.0) - exp(-sphere_radiance / settings.whitePoint.rbg * positional.exposure), vec3(1.0 / 2.2));

	color.rgb += diffuseColor.rgb * sphere_radiance;

    outColor.rgb = (color * ambientOcclusion) + emissiveColor.rgb;
    outColor.a = diffuseColor.a;
}
