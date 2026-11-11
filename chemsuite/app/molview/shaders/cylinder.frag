#version 440
// chemsuite 3D Viewer: ray-cylinder hit (open ends: atoms cover them), two-colour bonds
// (GPL-3.0-or-later).
layout(location = 0) in vec3 v_pos;
layout(location = 1) in vec4 v_start;
layout(location = 2) in vec4 v_axis;
layout(location = 3) in vec4 v_color0;
layout(location = 4) in vec4 v_color1;

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
    vec3 dir = normalize(v_pos);
    vec3 a = v_start.xyz;
    float r = v_start.w;
    vec3 ax = v_axis.xyz;
    float len = v_axis.w;
    vec3 oc = -a;                              // from the bond's start to the eye
    float dd = dot(dir, ax);
    float od = dot(oc, ax);
    vec3 dp = dir - ax * dd;
    vec3 op = oc - ax * od;
    float qa = dot(dp, dp);
    float qb = 2.0 * dot(dp, op);
    float qc = dot(op, op) - r * r;
    float disc = qb * qb - 4.0 * qa * qc;
    if (qa < 1e-10 || disc < 0.0)
        discard;
    float t = (-qb - sqrt(disc)) / (2.0 * qa);
    float s = od + t * dd;                     // how far along the bond the hit is
    if (s < 0.0 || s > len)
        discard;
    vec3 p = dir * t;
    vec3 n = normalize(p - (a + ax * s));
    vec4 clip = proj * vec4(p, 1.0);
    float z = clip.z / clip.w;
    gl_FragDepth = params.x > 0.5 ? z : z * 0.5 + 0.5;
    vec4 base = s < 0.5 * len ? v_color0 : v_color1;
    vec3 l = normalize(light.xyz);
    float diffuse = max(dot(n, l), 0.0);
    float spec = pow(max(dot(reflect(-l, n), -dir), 0.0), 48.0);
    vec3 shade = base.rgb * (0.32 + 0.68 * diffuse) + vec3(0.28 * spec);
    float fog = clamp((-p.z - params.y) / max(params.z - params.y, 1e-4), 0.0, 1.0);
    fragColor = vec4(mix(shade, background.rgb, fog * 0.65), base.a);
}
