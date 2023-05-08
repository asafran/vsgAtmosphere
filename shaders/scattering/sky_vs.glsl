#version 450
#extension GL_ARB_separate_shader_objects : enable

layout(push_constant) uniform PushConstants {
    mat4 projection;
    mat4 modelView;
} pc; //inverse

layout(location = 0) in vec2 osg_Vertex;
layout(location = 0) out vec3 outRay;

out gl_PerVertex{ vec4 gl_Position; };

void main()
{
    vec4 vertex = vec4(osg_Vertex, 0.0, 1.0);
    outRay = (pc.modelView * vec4((pc.projection * vertex).xyz, 0.0)).xyz;
    gl_Position = vertex;
}
