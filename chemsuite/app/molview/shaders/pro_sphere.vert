#version 440
// chemsuite protein view: one atom as a ray-traced sphere on a billboard (GPL-3.0-or-later).
// Geometry and look come from two buffers, so a selection uploads only the look.
layout(location = 0) in vec2 corner;      // per vertex: (-1,-1) .. (1,1)
layout(location = 1) in vec4 sphere;      // per instance: centre xyz, radius
layout(location = 2) in vec4 color;       // per instance: rgba
layout(location = 3) in float flags;      // per instance: selected 1 + hover 2 + ghost 4

layout(std140, binding = 0) uniform Frame {
    mat4 view;
    mat4 proj;
    vec4 light;
    vec4 params;      // x: 1 when clip depth runs 0..1; y: fog start; z: fog end; w: 1 orthographic
    vec4 background;
    vec4 slab;        // x, y: view depths kept; z: ghost fade
    vec4 highlight;   // the selection colour
};

layout(location = 0) out vec3 v_pos;
layout(location = 1) out vec4 v_sphere;
layout(location = 2) out vec4 v_color;
layout(location = 3) out float v_flags;

void main()
{
    vec3 c = (view * vec4(sphere.xyz, 1.0)).xyz;
    float r = sphere.w;
    float d = max(length(c), r * 1.0001);
    bool ortho = params.w > 0.5;               // parallel rays: the billboard faces the screen
    vec3 toward = ortho ? vec3(0.0, 0.0, 1.0) : -c / d;
    float front = d - r;
    float half_size = ortho ? r * 1.05 + 1e-4
                            : front * r / sqrt(max(d * d - r * r, 1e-8)) * 1.05 + 1e-4;
    vec3 helper = abs(toward.y) < 0.99 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 right = normalize(cross(toward, helper));
    vec3 up = cross(right, toward);
    vec3 p = c + toward * r + (corner.x * right + corner.y * up) * half_size;
    v_pos = p;
    v_sphere = vec4(c, r);
    v_color = color;
    v_flags = flags;
    gl_Position = proj * vec4(p, 1.0);
}
