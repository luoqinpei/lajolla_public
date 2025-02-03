#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    Spectrum base_color = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Spectrum F_m = base_color + (1 - base_color) * pow(1 - abs(dot(dir_out, h)), 5);
    Real aspect = sqrt(1 - anisotropic * 0.9);
    Real ax = std::max(0.0001, pow(roughness, 2) / aspect);
    Real ay = std::max(0.0001, pow(roughness, 2) * aspect);
    // Project the half vector to the local shading frame
    Vector3 h_local = to_local(frame, h);
    // Project the incoming direction to the local shading frame
    Vector3 wi_local = to_local(frame, dir_in);
    // Project the outgoing direction to the local shading frame
    Vector3 wo_local = to_local(frame, dir_out);
    Real D_m = 1 / (c_PI * ax * ay * pow((h_local.x * h_local.x )/ (ax * ax) + (h_local.y * h_local.y) / (ay * ay) + h_local.z * h_local.z, 2));
    Real Lambda_in = (sqrt(1 + (ax * ax * wi_local.x * wi_local.x + ay * ay * wi_local.y * wi_local.y) / (wi_local.z * wi_local.z)) - 1) / 2;
    Real Lambda_out = (sqrt(1 + (ax * ax * wo_local.x * wo_local.x + ay * ay * wo_local.y * wo_local.y) / (wo_local.z * wo_local.z)) - 1) / 2;
    Real G_m = 1 / ((1 + Lambda_in) * (1 + Lambda_out));
    Spectrum f_metal = F_m * D_m * G_m / (4 * abs(dot(frame.n, dir_in)));
    return f_metal;
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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

    // Half vector
    Vector3 h = normalize(dir_in + dir_out);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - anisotropic * 0.9);
    Real ax = std::max(0.0001, pow(roughness, 2) / aspect);
    Real ay = std::max(0.0001, pow(roughness, 2) * aspect);
    // Project the half vector to the local shading frame
    Vector3 h_local = to_local(frame, h);
    // Project the incoming direction to the local shading frame
    Vector3 wi_local = to_local(frame, dir_in);
    Real D_m = 1 / (c_PI * ax * ay * pow((h_local.x * h_local.x )/ (ax * ax) + (h_local.y * h_local.y) / (ay * ay) + h_local.z * h_local.z, 2));
    Real Lambda_in = (sqrt(1 + (ax * ax * wi_local.x * wi_local.x + ay * ay * wi_local.y * wi_local.y) / (wi_local.z * wi_local.z)) - 1) / 2;
    Real G_in = 1 / (1 + Lambda_in) ;

    Real weight = D_m * G_in / (4 * abs(dot(frame.n, dir_in)));
    return weight;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!        // Lambertian sampling
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Compute anisotropic parameters
    Real aspect = sqrt(1 - anisotropic * 0.9);
    Real ax = std::max(0.0001, pow(roughness, 2) / aspect);
    Real ay = std::max(0.0001, pow(roughness, 2) * aspect);

    // Sample half vector in local space
    Real phi = 2 * c_PI * rnd_param_uv.x;
    Real cos_theta = sqrt((1 - rnd_param_uv.y) / (1 + (ax * ax - 1) * rnd_param_uv.y));
    Real sin_theta = sqrt(1 - cos_theta * cos_theta);
    Vector3 h_local = Vector3{
        sin_theta * cos(phi),
        sin_theta * sin(phi),
        cos_theta
    };

    // Transform half vector to world space
    Vector3 h = to_world(frame, h_local);

    // Compute outgoing direction
    Vector3 dir_out = 2 * dot(dir_in, h) * h - dir_in;

    // Ensure dir_out is above surface
    if (dot(dir_out, frame.n) < 0) {
        return {};
    }

    Real pdf = pdf_sample_bsdf_op{  
    dir_in,      
    dir_out,
    vertex,
    texture_pool,
    dir}(bsdf);
    return BSDFSampleRecord{dir_out, Real(1.0), pdf};
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
