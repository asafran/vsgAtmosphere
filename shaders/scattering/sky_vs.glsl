#version 450
#extension GL_ARB_separate_shader_objects : enable

layout(push_constant) uniform PushConstants {
    mat4 projection;
    mat4 modelView;
} pc; //inverse

const vec2 vertexInput[4] = {{-1.0f, -1.0f},
                        {1.0f, -1.0f},
                        {-1.0f, 1.0f},
                        {1.0f, 1.0f}};

layout(location = 0) out vec3 outRay;
layout(location = 1) out vec2 outUV;

out gl_PerVertex{ vec4 gl_Position; };

void main()
{
    outUV = vertexInput[gl_VertexIndex];
    vec4 vertex = vec4(outUV, 0.0, 1.0);
    outRay = (pc.modelView * vec4((pc.projection * vertex).xyz, 0.0)).xyz;
    gl_Position = vertex;
}
