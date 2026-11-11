#version 440
// chemsuite protein view: lit triangle meshes, cartoons and surfaces (GPL-3.0-or-later).
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec4 color;       // look, per vertex: rgba
layout(location = 3) in float flags;      // selected 1 + hover 2 + ghost 4

layout(std140, binding = 0) uniform Frame {
    mat4 view;
    mat4 proj;
    vec4 light;
    vec4 params;
    vec4 background;
    vec4 slab;
    vec4 highlight;
};

layout(location = 0) out vec3 v_pos;
layout(location = 1) out vec3 v_normal;
layout(location = 2) out vec4 v_color;
layout(location = 3) out float v_flags;

void main()
{
    vec4 p = view * vec4(position, 1.0);
    v_pos = p.xyz;
    v_normal = mat3(view) * normal;
    v_color = color;
    v_flags = flags;
    gl_Position = proj * p;
}
