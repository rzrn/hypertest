in  vec4  _color;
in  vec3  _vertex;
out vec4  color;
out float fogFactor;

void main() {
    vec4 v1 = model(_vertex);
    vec4 v2 = view * v1;

    gl_Position = projection * v2;
    fogFactor   = max(getHorizFogFactor(length(v2.xyz / v2.w)), getVertFogFactor(v1.y / v1.w));
    color       = _color;
}
