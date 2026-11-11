#version 440
// chemsuite 3D Viewer: one atom as a ray-traced sphere on a billboard (GPL-3.0-or-later).
// The quad sits at the sphere's front, just large enough to cover its silhouette.
layout(location = 0) in vec2 corner;      // per vertex: (-1,-1) .. (1,1)
layout(location = 1) in vec4 sphere;      // per instance: centre xyz, radius
layout(location = 2) in vec4 color;       // per instance: rgba

layout(std140, binding = 0) uniform Frame {
    mat4 view;
    mat4 proj;
    vec4 light;       // view-space direction towards the light
    vec4 params;      // x: 1 when clip depth runs 0..1; y: fog start; z: fog end (view depth)
    vec4 background;  // fog colour
};

layout(location = 0) out vec3 v_pos;
layout(location = 1) out vec4 v_sphere;
layout(location = 2) out vec4 v_color;

void main()
{
    vec3 c = (view * vec4(sphere.xyz, 1.0)).xyz;
    float r = sphere.w;
    float d = max(length(c), r * 1.0001);
    vec3 toward = -c / d;
    float front = d - r;
    float half_size = front * r / sqrt(max(d * d - r * r, 1e-8)) * 1.05 + 1e-4;
    vec3 helper = abs(toward.y) < 0.99 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 right = normalize(cross(toward, helper));
    vec3 up = cross(right, toward);
    vec3 p = c + toward * r + (corner.x * right + corner.y * up) * half_size;
    v_pos = p;
    v_sphere = vec4(c, r);
    v_color = color;
    gl_Position = proj * vec4(p, 1.0);
}
