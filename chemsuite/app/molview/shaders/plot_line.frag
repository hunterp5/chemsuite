#version 440
// chemsuite 3D plots: lines in a flat colour (GPL-3.0-or-later).
layout(location = 0) in vec4 v_color;

layout(location = 0) out vec4 fragColor;

void main()
{
    fragColor = v_color;
}
