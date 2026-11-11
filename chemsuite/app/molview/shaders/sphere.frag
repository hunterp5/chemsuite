#version 440
// chemsuite 3D Viewer: ray-sphere hit, depth and shading (GPL-3.0-or-later).
layout(location = 0) in vec3 v_pos;
layout(location = 1) in vec4 v_sphere;
layout(location = 2) in vec4 v_color;

layout(std140, binding = 0) uniform Frame {
    mat4 view;
    mat4 proj;
    vec4 light;
    vec4 params;
    vec4 background;
};

layout(location = 0) out vec4 fragColor;

void main()
{
    vec3 dir = normalize(v_pos);              // the eye is at the origin of view space
    vec3 c = v_sphere.xyz;
    float r = v_sphere.w;
    float b = dot(dir, c);
    float h = b * b - dot(c, c) + r * r;
    if (h < 0.0)
        discard;
    vec3 p = dir * (b - sqrt(h));
    vec3 n = normalize(p - c);
    vec4 clip = proj * vec4(p, 1.0);
    float z = clip.z / clip.w;
    gl_FragDepth = params.x > 0.5 ? z : z * 0.5 + 0.5;
    vec3 l = normalize(light.xyz);
    float diffuse = max(dot(n, l), 0.0);
    float spec = pow(max(dot(reflect(-l, n), -dir), 0.0), 48.0);
    vec3 shade = v_color.rgb * (0.32 + 0.68 * diffuse) + vec3(0.28 * spec);
    float fog = clamp((-p.z - params.y) / max(params.z - params.y, 1e-4), 0.0, 1.0);
    fragColor = vec4(mix(shade, background.rgb, fog * 0.65), v_color.a);
}
