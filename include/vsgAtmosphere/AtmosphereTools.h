#ifndef ATMOSPHERETOOLS_H
#define ATMOSPHERETOOLS_H

#include "AtmoshpereConstants.h"
#include <cstdint>
#include <vector>
#include <assert.h>
#include <vsg/maths/mat4.h>
#include <vsg/maths/vec3.h>
#include <vsg/state/ImageInfo.h>
#include <vsg/state/Descriptor.h>
#include <vsg/state/DescriptorImage.h>

namespace atmosphere {
/*
vsg::ref_ptr<vsg::ImageInfo> createCubemap(uint32_t size)
{
    vsg::ref_ptr<vsg::Image> image = vsg::Image::create();
    image->usage |= (VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_STORAGE_BIT);
    image->format = VK_FORMAT_R32G32B32A32_SFLOAT;
    image->mipLevels = 1;
    image->extent = VkExtent3D{size, size, 1};
    image->imageType = VK_IMAGE_TYPE_2D;
    image->arrayLayers = 6;

    auto imageView = vsg::ImageView::create(image, VK_IMAGE_ASPECT_COLOR_BIT);
    imageView->viewType = VK_IMAGE_VIEW_TYPE_CUBE;
    auto sampler = vsg::Sampler::create();
    return vsg::ImageInfo::create(sampler, imageView, VK_IMAGE_LAYOUT_GENERAL);
}

void addImage(vsg::ref_ptr<vsg::ImageInfo> image, vsg::Descriptors &descriptors)
{
    int last = descriptors.empty() ? 0 : descriptors.back()->dstBinding + 1;
    auto storage = vsg::DescriptorImage::create(image, last++, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
    //auto sampled = vsg::DescriptorImage::create(image, last, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);
    descriptors.push_back(storage);
    //descriptors.push_back(sampled);
}

vsg::ref_ptr<vsg::ImageInfo> generateTexture(vsg::ref_ptr<vsg::Device> device, VkExtent3D extent, vsg::ref_ptr<vsg::Sampler> sampler = {}, bool init = false, VkFormat format = VK_FORMAT_R32G32B32A32_SFLOAT)
{
    vsg::ref_ptr<vsg::Image> image = vsg::Image::create();
    image->usage |= (VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT);
    image->format = format;
    image->mipLevels = 1;
    image->extent = extent;
    image->imageType = VK_IMAGE_TYPE_2D;
    image->arrayLayers = 1;

    if(extent.depth > 1)
    {
        image->imageType = VK_IMAGE_TYPE_3D;
        if(init)
            image->data = vsg::vec4Array3D::create(extent.width, extent.height, extent.depth, vsg::vec4{0.0f, 0.0f, 0.0f, 1.0f}, vsg::Data::Properties{format});
    }
    else
    {
        image->imageType = VK_IMAGE_TYPE_2D;
        if(init)
            image->data = vsg::vec4Array2D::create(extent.width, extent.height, vsg::vec4{0.0f, 0.0f, 0.0f, 1.0f}, vsg::Data::Properties{format});
    }

    image->compile(device);
    image->allocateAndBindMemory(device, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);

    auto imageView = vsg::ImageView::create(image, VK_IMAGE_ASPECT_COLOR_BIT);
    if(!sampler)
        sampler = vsg::Sampler::create();
    return vsg::ImageInfo::create(sampler, imageView, VK_IMAGE_LAYOUT_GENERAL);
}

vsg::ref_ptr<vsg::ImageInfo> generate3D(vsg::ref_ptr<vsg::Device> device, uint32_t width, uint32_t height, uint32_t depth, bool init = false)
{
    vsg::ref_ptr<vsg::Image> image = vsg::Image::create();
    image->usage |= (VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT);
    image->format = VK_FORMAT_R32G32B32A32_SFLOAT;
    image->mipLevels = 1;
    image->extent = VkExtent3D{width, height, depth};
    image->imageType = VK_IMAGE_TYPE_3D;
    image->arrayLayers = 1;

    if(init)


    image->compile(device);
    image->allocateAndBindMemory(device, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);

    auto imageView = vsg::ImageView::create(image, VK_IMAGE_ASPECT_COLOR_BIT);
    auto sampler = vsg::Sampler::create();
    return vsg::ImageInfo::create(sampler, imageView, VK_IMAGE_LAYOUT_GENERAL);
}

vsg::ref_ptr<vsg::ImageInfo> generateNoise(vsg::ref_ptr<vsg::Device> device, uint32_t size)
{
    vsg::ref_ptr<vsg::Image> image = vsg::Image::create();
    image->usage |= (VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT);
    image->format = VK_FORMAT_R16G16B16A16_SFLOAT;
    image->mipLevels = 1;
    image->extent = VkExtent3D{size, size, size};
    image->imageType = VK_IMAGE_TYPE_3D;
    image->arrayLayers = 1;

    image->compile(device);
    image->allocateAndBindMemory(device, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);

    auto imageView = vsg::ImageView::create(image, VK_IMAGE_ASPECT_COLOR_BIT);
    auto sampler = vsg::Sampler::create();
    return vsg::ImageInfo::create(sampler, imageView, VK_IMAGE_LAYOUT_GENERAL);
}

vsg::ref_ptr<vsg::Data> convertArray(vsg::ref_ptr<vsg::Data> array)
{
    auto width = std::sqrt(array->width());
    auto height = array->height();
    auto data = static_cast<vsg::ubvec4*>(array->dataRelease());
    if(!data || !array->is_compatible(typeid(vsg::ubvec4Array2D)))
        return {};
    return vsg::ubvec4Array3D::create(height, width, width, data, array->properties);
}
*/
inline double interpolate(
        const std::vector<double>& wavelengths,
        const std::vector<double>& wavelength_function,
        double wavelength)
{
    assert(wavelength_function.size() == wavelengths.size());
    if (wavelength < wavelengths[0]) {
      return wavelength_function[0];
    }
    for (unsigned int i = 0; i < wavelengths.size() - 1; ++i) {
      if (wavelength < wavelengths[i + 1]) {
        double u =
            (wavelength - wavelengths[i]) / (wavelengths[i + 1] - wavelengths[i]);
        return
            wavelength_function[i] * (1.0 - u) + wavelength_function[i + 1] * u;
      }
    }
    return wavelength_function[wavelength_function.size() - 1];
}

constexpr double cieColorMatchingFunctionTableValue(double wavelength, int column)
{
    if (wavelength <= kLambdaMin || wavelength >= kLambdaMax)
        return 0.0;

    double u = (wavelength - kLambdaMin) / 5.0;
    int row = static_cast<int>(std::floor(u));
    assert(row >= 0 && row + 1 < 95);
    assert(CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * row] <= wavelength &&
           CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * (row + 1)] >= wavelength);
    u -= row;
    return CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * row + column] * (1.0 - u) +
        CIE_2_DEG_COLOR_MATCHING_FUNCTIONS[4 * (row + 1) + column] * u;
}

inline vsg::vec4 computeSpectralRadianceToLuminanceFactors(
    const std::vector<double>& wavelengths,
    const std::vector<double>& solar_irradiance,
    double lambda_power) {

  auto k_r = 0.0;
  auto k_g = 0.0;
  auto k_b = 0.0;

  double solar_r = interpolate(wavelengths, solar_irradiance, kLambdaR);
  double solar_g = interpolate(wavelengths, solar_irradiance, kLambdaG);
  double solar_b = interpolate(wavelengths, solar_irradiance, kLambdaB);
  int dlambda = 1;
  for (int lambda = kLambdaMin; lambda < kLambdaMax; lambda += dlambda) {
      double x_bar = cieColorMatchingFunctionTableValue(lambda, 1);
      double y_bar = cieColorMatchingFunctionTableValue(lambda, 2);
      double z_bar = cieColorMatchingFunctionTableValue(lambda, 3);
    const double* xyz2srgb = XYZ_TO_SRGB;
    double r_bar =
        xyz2srgb[0] * x_bar + xyz2srgb[1] * y_bar + xyz2srgb[2] * z_bar;
    double g_bar =
        xyz2srgb[3] * x_bar + xyz2srgb[4] * y_bar + xyz2srgb[5] * z_bar;
    double b_bar =
        xyz2srgb[6] * x_bar + xyz2srgb[7] * y_bar + xyz2srgb[8] * z_bar;
    double irradiance = interpolate(wavelengths, solar_irradiance, lambda);
    k_r += r_bar * irradiance / solar_r *
        pow(lambda / kLambdaR, lambda_power);
    k_g += g_bar * irradiance / solar_g *
        pow(lambda / kLambdaG, lambda_power);
    k_b += b_bar * irradiance / solar_b *
        pow(lambda / kLambdaB, lambda_power);
  }
  k_r *= MAX_LUMINOUS_EFFICACY * dlambda;
  k_g *= MAX_LUMINOUS_EFFICACY * dlambda;
  k_b *= MAX_LUMINOUS_EFFICACY * dlambda;

  return {static_cast<float>(k_r), static_cast<float>(k_g), static_cast<float>(k_b), 0.0f};
}

inline vsg::vec4 toVector(const std::vector<double>& wavelengths, const std::vector<double>& v, const vsg::vec3& lambdas, double scale)
{
    auto r = static_cast<float>(interpolate(wavelengths, v, lambdas[0]) * scale);
    auto g = static_cast<float>(interpolate(wavelengths, v, lambdas[1]) * scale);
    auto b = static_cast<float>(interpolate(wavelengths, v, lambdas[2]) * scale);
    return {r, g, b, 1.0f};
}

constexpr vsg::mat4 toMatrix(double arr[])
{
    return vsg::mat4(
        (float)arr[0], (float)arr[3], (float)arr[6], 0.0f,
        (float)arr[1], (float)arr[4], (float)arr[7], 0.0f,
        (float)arr[2], (float)arr[5], (float)arr[8], 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f);
}

constexpr double coeff(double dlambda, double lambda, int component)
{
    // Note that we don't include MAX_LUMINOUS_EFFICACY here, to avoid
    // artefacts due to too large values when using half precision on GPU.
    // We add this term back in kAtmosphereShader, via
    // SKY_SPECTRAL_RADIANCE_TO_LUMINANCE (see also the comments in the
    // Model constructor).
    double x = cieColorMatchingFunctionTableValue(lambda, 1);
    double y = cieColorMatchingFunctionTableValue(lambda, 2);
    double z = cieColorMatchingFunctionTableValue(lambda, 3);
    return static_cast<float>((
            XYZ_TO_SRGB[component * 3] * x +
            XYZ_TO_SRGB[component * 3 + 1] * y +
            XYZ_TO_SRGB[component * 3 + 2] * z) * dlambda);
}

constexpr double day2000(tm time) {
    int d1 = 0;
    int b = 0;
    int c = 0;
    int greg = 0;

    auto y = time.tm_year;
    auto m = time.tm_mon;
    auto d = time.tm_mday;
    auto h = time.tm_hour;

    greg = y * 10000 + m * 100 + d;
    if (m == 1 || m == 2) {
        y = y - 1;
        m = m + 12;
    }
    //  reverts to Julian calendar before 4th Oct 1582
    //  no good for UK, America or Sweeden!

    if (greg > 15821004) {
        auto a = std::floor(static_cast<double>(y) / 100.0);
        b = 2 - a + static_cast<int>(std::floor(a / 4.0));
    }
    else {
        b = 0;
    }
    c = static_cast<int>(std::floor(365.25 * y));
    d1 = static_cast<int>(std::floor(30.6001 * (m + 1)));
    return (b + c + d1 - 730550.5 + d + h / 24);
}
}

#endif // ATMOSPHERETOOLS_H
