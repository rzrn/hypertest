in  vec3  _vertex;
out float fogFactor;

void main() {
    vec4 v1 = model(_vertex);
    vec4 v2 = view * v1;

    gl_Position = projection * v2;
    fogFactor   = getVertFogFactor(v1.y / v1.w);
}
