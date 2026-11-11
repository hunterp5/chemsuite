#version 440
// chemsuite protein view: one bond as a ray-traced cylinder in a box (GPL-3.0-or-later).
layout(location = 0) in vec3 corner;      // per vertex: x 0..1 along the bond, y and z -1..1
layout(location = 1) in vec4 start;       // per instance: first end xyz, radius
layout(location = 2) in vec4 end;         // per instance: second end xyz, unused
layout(location = 3) in vec4 color0;      // look: the first half's colour
layout(location = 4) in vec4 color1;      // the second half's colour
layout(location = 5) in vec2 flags;       // each half's flags

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
layout(location = 1) out vec4 v_start;
layout(location = 2) out vec4 v_axis;
layout(location = 3) out vec4 v_color0;
layout(location = 4) out vec4 v_color1;
layout(location = 5) out vec2 v_flags;

void main()
{
    vec3 a = (view * vec4(start.xyz, 1.0)).xyz;
    vec3 b = (view * vec4(end.xyz, 1.0)).xyz;
    float r = start.w;
    vec3 axis = b - a;
    float len = max(length(axis), 1e-6);
    vec3 ax = axis / len;
    vec3 helper = abs(ax.y) < 0.99 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 u = normalize(cross(ax, helper));
    vec3 w = cross(ax, u);
    vec3 p = a + ax * (corner.x * len) + (u * corner.y + w * corner.z) * (r * 1.02);
    v_pos = p;
    v_start = vec4(a, r);
    v_axis = vec4(ax, len);
    v_color0 = color0;
    v_color1 = color1;
    v_flags = flags;
    gl_Position = proj * vec4(p, 1.0);
}
