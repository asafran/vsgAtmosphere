#version 450
#extension GL_ARB_separate_shader_objects : enable

layout(binding = 0, set = 0) uniform samplerCube envMap;
layout(location = 0) in vec3 inUVW;
layout(location = 0) out vec4 outColor;

void main()
{
    outColor = textureLod(envMap, inUVW, 0);
}
