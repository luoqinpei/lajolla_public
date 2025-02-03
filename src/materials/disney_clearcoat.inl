#include "../microfacet.h"
#include <random>

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Half vector
    Vector3 h = normalize(dir_in + dir_out);
    // Project the half vector to the local shading frame
    Vector3 h_local = to_local(frame, h);
    Real eta = 1.5;
    Real R0 = pow((eta - 1) / (eta + 1), 2);
    Real F_c = R0 + (1 - R0) * pow(1 - abs(dot(h, dir_out)), 5);
    Real alpha_g = (1 - clearcoat) * 0.1 + clearcoat * 0.001;
    Real D_c = (alpha_g * alpha_g - 1) / (c_PI * log(alpha_g * alpha_g) *(1 + (alpha_g * alpha_g - 1) * (h_local.z * h_local.z)));

    // Project the incoming direction to the local shading frame
    Vector3 wi_local = to_local(frame, dir_in);
    // Project the outgoing direction to the local shading frame
    Vector3 wo_local = to_local(frame, dir_out);
    Real Lambda_in = (sqrt(1 + wi_local.x * wi_local.x * 0.25 * 0.25 / (wi_local.z * wi_local.z) 
                        + wi_local.y * wi_local.y * 0.25 * 0.25 / (wi_local.z * wi_local.z)) - 1) / 2;
    Real Lambda_out = (sqrt(1 + wo_local.x * wo_local.x * 0.25 * 0.25 / (wo_local.z * wo_local.z) 
                        + wo_local.y * wo_local.y * 0.25 * 0.25 / (wo_local.z * wo_local.z)) - 1) / 2;
    Real G_c = 1 / ((1 + Lambda_in) * (1 + Lambda_out));
    Spectrum f_clearcoat;
    f_clearcoat.x = F_c * D_c * G_c / (4 * abs(dot(frame.n, dir_in)));
    f_clearcoat.y = F_c * D_c * G_c / (4 * abs(dot(frame.n, dir_in)));
    f_clearcoat.z = F_c * D_c * G_c / (4 * abs(dot(frame.n, dir_in)));
    return f_clearcoat;
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    // Importance sampling the clearcoat lobe by choosing a micro normal proportional to D_c

    // Genearte 2D random number 
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(0.0, 1.0);

    Real u0 = dis(gen);
    Real u1 = dis(gen);
    Real clearcoat = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real alpha_g = (1 - clearcoat) * 0.1 + clearcoat * 0.001;
    Real cos_el = sqrt((1 - pow(alpha_g * alpha_g, 1 - u0))/ (1 - alpha_g * alpha_g));
    Real sin_el = sqrt(1 - cos_el * cos_el);
    Real h_az = 2 * c_PI * u1;
    Real h_x = cos(h_az) * sin_el;
    Real h_y = sin(h_az) * sin_el;
    Real h_z = cos_el;
    Vector3 h = Vector3(h_x, h_y, h_z);

    Real D_c = (alpha_g * alpha_g) / (c_PI * log(alpha_g * alpha_g) *(1 + (alpha_g * alpha_g - 1) * (h_z * h_z)));

    Real pdf = D_c * abs(dot(frame.n, h)) / (4 * abs(dot(h, dir_out)));
    return pdf;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
   // Genearte 2D random number 
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(0.0, 1.0);

    Real u0 = dis(gen);
    Real u1 = dis(gen);
    Real clearcoat = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real alpha_g = (1 - clearcoat) * 0.1 + clearcoat * 0.001;
    Real cos_el = sqrt((1 - pow(alpha_g * alpha_g, 1 - u0))/ (1 - alpha_g * alpha_g));
    Real sin_el = sqrt(1 - cos_el * cos_el);
    Real h_az = 2 * c_PI * u1;
    Real h_x = cos(h_az) * sin_el;
    Real h_y = sin(h_az) * sin_el;
    Real h_z = cos_el;
    Vector3 h = Vector3(h_x, h_y, h_z);

    // Compute outgoing direction
    Vector3 dir_out = 2 * dot(dir_in, h) * h - dir_in;
    Real D_c = (alpha_g * alpha_g) / (c_PI * log(alpha_g * alpha_g) *(1 + (alpha_g * alpha_g - 1) * (h_z * h_z)));

    Real pdf = D_c * abs(dot(frame.n, h)) / (4 * abs(dot(h, dir_out)));

    return BSDFSampleRecord{dir_out, Real(1.0), pdf};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
