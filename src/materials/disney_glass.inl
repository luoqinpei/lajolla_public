#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // Half vector
    Spectrum base_color = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Vector3 h, h_local;
    // Project the incoming direction to the local shading frame
    Vector3 wi_local = to_local(frame, dir_in);
    // Project the outgoing direction to the local shading frame
    Vector3 wo_local = to_local(frame, dir_out);

    if (reflect) {
        h = normalize(dir_in + dir_out);
        h_local = normalize(wi_local + wo_local);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
        h_local = normalize(wi_local + wo_local * eta);
    }
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    Real h_dot_in = dot(h, dir_in);
    Real h_dot_out = dot(h, dir_out);

    Real R_s = (h_dot_in - eta * h_dot_out) / (h_dot_in + eta * h_dot_out);
    Real R_p = (eta * h_dot_in - h_dot_out) / (eta * h_dot_in + h_dot_out);
    Real F_g = (R_s * R_s + R_p * R_p) / 2;
    // Real F_g = fresnel_dielectric(h_dot_in, eta);

    Real aspect = sqrt(1 - anisotropic * 0.9);
    Real ax = std::max(0.0001, pow(roughness, 2) / aspect);
    Real ay = std::max(0.0001, pow(roughness, 2) * aspect);
    Real D_g = 1 / (c_PI * ax * ay * pow((h_local.x * h_local.x )/ (ax * ax) + (h_local.y * h_local.y) / (ay * ay) + h_local.z * h_local.z, 2));
    Real Lambda_in = (sqrt(1 + (ax * ax * wi_local.x * wi_local.x + ay * ay * wi_local.y * wi_local.y) / (wi_local.z * wi_local.z)) - 1) / 2;
    Real Lambda_out = (sqrt(1 + (ax * ax * wo_local.x * wo_local.x + ay * ay * wo_local.y * wo_local.y) / (wo_local.z * wo_local.z)) - 1) / 2;
    Real G_g = 1 / ((1 + Lambda_in) * (1 + Lambda_out));
    Spectrum f_glass;
    if (reflect) {
        f_glass = base_color * F_g * D_g * G_g / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        f_glass = sqrt(base_color) * (1 - F_g) * D_g * G_g * fabs(dot(h, dir_out) * dot(h, dir_in)) / (fabs(dot(frame.n, dir_in)) * pow(dot(h, dir_in) + eta * dot(h, dir_out), 2));
    }
    return f_glass;
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    
     // Half vector
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);
    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
    }
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real R_s = (dot(h, dir_in) - eta * dot(h, dir_out)) / (dot(h, dir_in) + eta * dot(h, dir_out));
    Real R_p = (eta * dot(h, dir_in) - dot(h, dir_out)) / (eta * dot(h, dir_in) + dot(h, dir_out));
    Real F_g = (R_s * R_s + R_p * R_p) / 2;

    // Real F_g = fresnel_dielectric(dot(h, dir_in), eta);

    Real aspect = sqrt(1 - anisotropic * 0.9);
    Real ax = std::max(0.0001, pow(roughness, 2) / aspect);
    Real ay = std::max(0.0001, pow(roughness, 2) * aspect);
    // Project the half vector to the local shading frame
    Vector3 h_local = to_local(frame, h);
    // Project the incoming direction to the local shading frame
    Vector3 wi_local = to_local(frame, dir_in);
    Real D_g = 1 / (c_PI * ax * ay * pow((h_local.x * h_local.x )/ (ax * ax) + (h_local.y * h_local.y) / (ay * ay) + h_local.z * h_local.z, 2));
    Real Lambda_in = (sqrt(1 + (ax * ax * wi_local.x * wi_local.x + ay * ay * wi_local.y * wi_local.y) / (wi_local.z * wi_local.z)) - 1) / 2;
    Real G_in = 1 / (1 + Lambda_in);   

    if (reflect) {
        return (F_g * D_g * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        return (1 - F_g) * D_g * G_in * fabs(dot(h, dir_out) * dot(h, dir_in)) / (fabs(dot(frame.n, dir_in)) * pow(dot(h, dir_in) + eta * dot(h, dir_out), 2));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect  = std::sqrt(1.0 - 0.9 * anisotropic);
    Real alpha_x = std::max(Real(0.0001), (roughness * roughness) / aspect);
    Real alpha_y = std::max(Real(0.0001), (roughness * roughness) * aspect);

    Vector3 local_dir_in = to_local(frame, dir_in);
    if (local_dir_in.z < 0) {
        local_dir_in.z = -local_dir_in.z;
    }
    Vector3 v_hemi = normalize(Vector3(alpha_x * local_dir_in.x, alpha_y * local_dir_in.y, local_dir_in.z));

    Real r = std::sqrt(rnd_param_uv.x);
    Real phi = 2 * c_PI * rnd_param_uv.y;
    Real t1 = r * std::cos(phi);
    Real t2 = r * std::sin(phi);
    Real s = 0.5 * (1.0 + v_hemi.z);
    t2 = (1.0 - s) * std::sqrt(std::max(0.0, 1.0 - t1 * t1)) + s * t2;
    Vector3 diskN(t1, t2, std::sqrt(std::max(0.0, 1.0 - t1 * t1 - t2 * t2)));
    Frame hemiFrame(v_hemi);
    Vector3 n_hemi = to_world(hemiFrame, diskN);

    Vector3 local_micro_normal = normalize(Vector3(alpha_x * n_hemi.x, alpha_y * n_hemi.y, std::max(0.0, n_hemi.z)));

    Vector3 h = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h= -h;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(h, dir_in);
    Real F_g = fresnel_dielectric(dot(h, dir_in), eta);

    if (rnd_param_w <= F_g) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (h_dot_in < 0) {
            h = -h;
        }
        Real h_dot_out= sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * h;
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
