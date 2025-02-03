Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!
    // Half vector
    Vector3 h = normalize(dir_in + dir_out);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1)); 
    Real subsurface = eval(
        bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real F_D90 = 0.5 + 2 * roughness * pow(abs(dot(h, dir_out)), 2);
    Real F_D_in = 1 + (F_D90 - 1) * pow(1 - abs(dot(dir_in, frame.n)), 5);
    Real F_D_out = 1 + (F_D90 - 1) * pow(1 - abs(dot(dir_out, frame.n)), 5);
    Spectrum f_baseDiffuse = base_color / c_PI * F_D_in * F_D_out * abs(dot(frame.n, dir_out));
    Real F_SS90 = roughness * pow(abs(dot(h, dir_out)), 2);
    Real F_SS_in =  1 + (F_SS90 - 1) * pow(1 - abs(dot(dir_in, frame.n)), 5);
    Real F_SS_out = 1 + (F_SS90 - 1) * pow(1 - abs(dot(dir_out, frame.n)), 5);
    Spectrum f_subsurface = 1.25 * base_color / c_PI * (F_SS_in * F_SS_out * (1 / (abs(dot(dir_in, frame.n)) + abs(dot(dir_out, frame.n))) - 0.5) + 0.5) * abs(dot(frame.n, dir_out));

    Spectrum f_diffuse = (1 - subsurface) * f_baseDiffuse + subsurface * f_subsurface;
    
    return f_diffuse;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!

    // For Lambertian, we importance sample the cosine hemisphere domain.
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    
    // Homework 1: implement this!
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
