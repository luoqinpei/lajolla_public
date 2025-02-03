#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Spectrum base_color = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(
        bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(
        bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(
        bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(
        bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(
        bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = bsdf.eta;

    // Half vector
    Vector3 h = normalize(dir_in + dir_out);
    Vector3 h_local = to_local(frame, h);
    // F_m
    Spectrum C_tint;
    if (luminance(base_color) > 0) {
        C_tint = base_color/luminance(base_color);
    }
    else {
        C_tint = make_const_spectrum(1);
    }
    Spectrum K_s = (1 - specular_tint) + specular_tint * C_tint;
    Real R0 = pow((eta - 1) / (eta + 1), 2);
    Spectrum C_0 = specular * R0 * (1 - metallic) * K_s + metallic * base_color;
    Spectrum F_m_hat = reflect > 0 ? C_0 + (1 - C_0) * pow(1 - dot(h, dir_out), 5) : make_zero_spectrum();

    // Diffuse   // Genearte 2D random number 
    Real F_D90 = 0.5 + 2 * roughness * pow(abs(dot(h, dir_out)), 2);
    Real F_D_in = 1 + (F_D90 - 1) * pow(1 - abs(dot(dir_in, frame.n)), 5);
    Real F_D_out = 1 + (F_D90 - 1) * pow(1 - abs(dot(dir_out, frame.n)), 5);
    Spectrum f_baseDiffuse = base_color / c_PI * F_D_in * F_D_out * abs(dot(frame.n, dir_out));
    Real F_SS90 = roughness * pow(abs(dot(h, dir_out)), 2);
    Real F_SS_in =  1 + (F_SS90 - 1) * pow(1 - abs(dot(dir_in, frame.n)), 5);
    Real F_SS_out = 1 + (F_SS90 - 1) * pow(1 - abs(dot(dir_out, frame.n)), 5);
    Spectrum f_subsurface = 1.25 * base_color / c_PI * (F_SS_in * F_SS_out * (1 / (abs(dot(dir_in, frame.n)) + abs(dot(dir_out, frame.n))) - 0.5) + 0.5) * abs(dot(frame.n, dir_out));

    Spectrum f_diffuse = reflect > 0 ? (1 - subsurface) * f_baseDiffuse + subsurface * f_subsurface : make_zero_spectrum();
    
    // Sheen
    Spectrum C_sheen = (1 - sheen_tint) + sheen_tint * C_tint;
    Spectrum f_sheen = reflect > 0 ? C_sheen * pow(1 - fabs(dot(h, dir_out)), 5) * fabs(dot(frame.n, dir_out)) : make_zero_spectrum();

    // Clearcoat
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
    if (reflect <= 0)
        f_clearcoat = make_zero_spectrum();

    // glass
    Real h_dot_in = dot(h, dir_in);
    Real h_dot_out = dot(h, dir_out);

    Real R_s = (h_dot_in - eta * h_dot_out) / (h_dot_in + eta * h_dot_out);
    Real R_p = (eta * h_dot_in - h_dot_out) / (eta * h_dot_in + h_dot_out);
    Real F_g = (R_s * R_s + R_p * R_p) / 2;

    Real aspect = sqrt(1 - anisotropic * 0.9);
    Real ax = std::max(0.0001, pow(roughness, 2) / aspect);
    Real ay = std::max(0.0001, pow(roughness, 2) * aspect);
    Real D_g = 1 / (c_PI * ax * ay * pow((h_local.x * h_local.x )/ (ax * ax) + (h_local.y * h_local.y) / (ay * ay) + h_local.z * h_local.z, 2));
    Real Lambda_in_g = (sqrt(1 + (ax * ax * wi_local.x * wi_local.x + ay * ay * wi_local.y * wi_local.y) / (wi_local.z * wi_local.z)) - 1) / 2;
    Real Lambda_out_g = (sqrt(1 + (ax * ax * wo_local.x * wo_local.x + ay * ay * wo_local.y * wo_local.y) / (wo_local.z * wo_local.z)) - 1) / 2;
    Real G_g = 1 / ((1 + Lambda_in_g) * (1 + Lambda_out_g));
    Spectrum f_glass;
    if (reflect) {
        f_glass = base_color * F_g * D_g * G_g / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        f_glass = sqrt(base_color) * (1 - F_g) * D_g * G_g * fabs(dot(h, dir_out) * dot(h, dir_in)) / (fabs(dot(frame.n, dir_in)) * pow(dot(h, dir_in) + eta * dot(h, dir_out), 2));
    }
    // Metallic
    Real D_m = 1 / (c_PI * ax * ay * pow((h_local.x * h_local.x )/ (ax * ax) + (h_local.y * h_local.y) / (ay * ay) + h_local.z * h_local.z, 2));
    Real G_m = 1 / ((1 + Lambda_in_g) * (1 + Lambda_out_g));
    Spectrum f_metal_hat = F_m_hat * D_m * G_m / (4 * abs(dot(frame.n, dir_in)));

    Spectrum f_disney = (1 - specular_transmission) * (1 - metallic) * f_diffuse +
            (1 - metallic) * sheen * f_sheen +
            (1 - specular_transmission * (1 - metallic)) * f_metal_hat +
            0.25 * clearcoat * f_clearcoat +
            (1 - metallic) * specular_transmission * f_glass;

    return f_disney;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // Define the weight for each lobe
    Real specular_transmission = eval(
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(
        bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(
        bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(
        bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(
        bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(
        bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = bsdf.eta;

    Real diffuseWeight = (1 - specular_transmission) * (1 - metallic);
    Real metalWeight = 1 - specular_transmission * (1 - metallic);
    Real glassWeight = (1 - metallic) * specular_transmission;
    Real clearcoatWeight = 0.25 * clearcoat;

    Real norm = 1 / (diffuseWeight + metalWeight + glassWeight + clearcoatWeight);
    Real p_diffuse = norm * diffuseWeight;
    Real p_metal = norm * metalWeight;
    Real p_glass = norm * glassWeight;
    Real p_clearcoat = norm * clearcoatWeight;

    Real pdfs[4] = {0, 0, 0, 0};

    if (dot(vertex.geometric_normal, dir_in) <= 0) // light from inside the object
    {
        // only consider the glass lobes
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
    else
    {
        // Define the diffuse lobe
        pdfs[0] = fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
        // Define the metal lobe
        // Half vector
        Vector3 h = normalize(dir_in + dir_out);

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
        pdfs[1] = D_m * G_in / (4 * abs(dot(frame.n, dir_in)));

        // Glass lobe
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

        Real R_s = (dot(h, dir_in) - eta * dot(h, dir_out)) / (dot(h, dir_in) + eta * dot(h, dir_out));
        Real R_p = (eta * dot(h, dir_in) - dot(h, dir_out)) / (eta * dot(h, dir_in) + dot(h, dir_out));
        Real F_g = (R_s * R_s + R_p * R_p) / 2;

        Real D_g = 1 / (c_PI * ax * ay * pow((h_local.x * h_local.x )/ (ax * ax) + (h_local.y * h_local.y) / (ay * ay) + h_local.z * h_local.z, 2));

        if (reflect) {
            pdfs[2] = (F_g * D_g * G_in) / (4 * fabs(dot(frame.n, dir_in)));
        } else {
            pdfs[2] = (1 - F_g) * D_g * G_in * fabs(dot(h, dir_out) * dot(h, dir_in)) / (fabs(dot(frame.n, dir_in)) * pow(dot(h, dir_in) + eta * dot(h, dir_out), 2));
        }

        // Clearcoat lobe
        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_real_distribution<> dis(0.0, 1.0);

        Real u0 = dis(gen);
        Real u1 = dis(gen);

        Real alpha_g = (1 - clearcoat) * 0.1 + clearcoat * 0.001;
        Real cos_el = sqrt((1 - pow(alpha_g * alpha_g, 1 - u0))/ (1 - alpha_g * alpha_g));
        Real sin_el = sqrt(1 - cos_el * cos_el);
        Real h_az = 2 * c_PI * u1;
        Real h_x = cos(h_az) * sin_el;
        Real h_y = sin(h_az) * sin_el;
        Real h_z = cos_el;
        Vector3 hh = Vector3(h_x, h_y, h_z);

        Real D_c = (alpha_g * alpha_g) / (c_PI * log(alpha_g * alpha_g) *(1 + (alpha_g * alpha_g - 1) * (h_z * h_z)));

        pdfs[3] = D_c * abs(dot(frame.n, hh)) / (4 * abs(dot(hh, dir_out)));

        return pdfs[0] * p_diffuse + pdfs[1] * p_metal + pdfs[2] * p_glass + pdfs[3] * p_clearcoat;
    }

}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Define the weight for each lobe
    Real specular_transmission = eval(
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(
        bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(
        bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(
        bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(
        bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(
        bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = bsdf.eta;

    Real diffuseWeight = (1 - specular_transmission) * (1 - metallic);
    Real metalWeight = 1 - specular_transmission* (1 - metallic);
    Real glassWeight = (1 - metallic) * specular_transmission;
    Real clearcoatWeight = 0.25 * clearcoat;

    // Compute anisotropic parameters
    Real aspect = sqrt(1 - anisotropic * 0.9);
    Real ax = std::max(0.0001, pow(roughness, 2) / aspect);
    Real ay = std::max(0.0001, pow(roughness, 2) * aspect);

    if (dot(vertex.geometric_normal, dir_in) <= 0) // light from inside the object
    {
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

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
    else
    {
        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_real_distribution<> dis(0.0, 1.0);
        // Sample a lobe based on the weights
        // Randomly choose a lobe
        Real sum = diffuseWeight + metalWeight + glassWeight + clearcoatWeight;
        Real u = dis(gen) * sum;
        if (u < diffuseWeight) {
            return BSDFSampleRecord{
                to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
                Real(0) /* eta */, Real(1) /* roughness */};
        } else if (u < diffuseWeight + metalWeight) {
            
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

            // Compute PDF and record
            Real D = 1 / (c_PI * ax * ay * pow(
                (h_local.x * h_local.x) / (ax * ax) +
                (h_local.y * h_local.y) / (ay * ay) +
                h_local.z * h_local.z, 2));
            Real pdf = D * abs(h_local.z) / (4 * abs(dot(dir_in, h)));

            return BSDFSampleRecord{dir_out, Real(1.0), pdf};
            
        } else if (u < diffuseWeight + metalWeight + glassWeight) {

            Vector3 local_dir_in = to_local(frame, dir_in);
            if (local_dir_in.z < 0) {
                local_dir_in.z = -local_dir_in.z;
            }
            Vector3 v_hemi = normalize(Vector3(ax * local_dir_in.x, ay * local_dir_in.y, local_dir_in.z));

            Real r = std::sqrt(rnd_param_uv.x);
            Real phi = 2 * c_PI * rnd_param_uv.y;
            Real t1 = r * std::cos(phi);
            Real t2 = r * std::sin(phi);
            Real s = 0.5 * (1.0 + v_hemi.z);
            t2 = (1.0 - s) * std::sqrt(std::max(0.0, 1.0 - t1 * t1)) + s * t2;
            Vector3 diskN(t1, t2, std::sqrt(std::max(0.0, 1.0 - t1 * t1 - t2 * t2)));
            Frame hemiFrame(v_hemi);
            Vector3 n_hemi = to_world(hemiFrame, diskN);

            Vector3 local_micro_normal = normalize(Vector3(ax * n_hemi.x, ay * n_hemi.y, std::max(0.0, n_hemi.z)));

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

        } else {
        // Genearte 2D random number 
            std::random_device rd;
            std::mt19937 gen(rd());

            std::uniform_real_distribution<> dis(0.0, 1.0);

            Real u0 = dis(gen);
            Real u1 = dis(gen);

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
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
