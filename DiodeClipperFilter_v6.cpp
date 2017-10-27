/**
 * \file DiodeClipperFilter_v6.cpp
 */

#include "DiodeClipperFilter_v6.h"

#include <boost/math/special_functions/sign.hpp>

#include <ATK/EQ/ButterworthFilter.h>
#include <ATK/EQ/IIRFilter.h>

#include <ATK/Distortion/DiodeClipperFilter.h>

#include <ATK/Mock/SimpleSinusGeneratorFilter.h>

#include <ATK/Tools/DecimationFilter.h>
#include <ATK/Tools/OversamplingFilter.h>

#include <ATK/Utility/fmath.h>
#include <ATK/Utility/ScalarNewtonRaphson.h>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#define PROCESSSIZE (1024)

namespace ATKADC
{
  template<typename DataType_>
  class DiodeClipperFilterv6<DataType_>::SimpleOverdriveFunction
  {
  public:
    typedef DataType_ DataType;
  protected:
    const DataType R;
    const DataType C;
    DataType is;
    DataType vt;

    DataType ieq;

    DataType expdiode_y1_p;
    DataType expdiode_y1_m;
  public:
    SimpleOverdriveFunction(DataType dt, DataType R, DataType C, DataType is, DataType vt)
    :R(R), C(2 * C / dt), is(is), vt(vt), ieq(0), expdiode_y1_p(1), expdiode_y1_m(1)
    {
    }
    
    std::pair<DataType, DataType> operator()(const DataType* ATK_RESTRICT input, DataType* ATK_RESTRICT output, DataType y1)
    {
      auto x1 = input[0];
      expdiode_y1_p = fmath::exp(y1 / vt);
      expdiode_y1_m = 1 / expdiode_y1_p;
	        
      std::pair<DataType, DataType> diode = std::make_pair(is * (expdiode_y1_p - expdiode_y1_m), is * (expdiode_y1_p + expdiode_y1_m) / vt);
      return std::make_pair(y1 * (1 + R * C) - x1 + R * diode.first - ieq * R, (1 + R * C) + R * diode.second);
    }

    void update_state(const DataType* ATK_RESTRICT input, DataType* ATK_RESTRICT output)
    {
      auto y1 = output[0];
      ieq = 2 * C * y1 - ieq;
    }

    DataType estimate(const DataType* ATK_RESTRICT input, DataType* ATK_RESTRICT output)
    {
      auto x0 = input[-1];
      auto x1 = input[0];
      auto y0 = output[-1];
      return affine_estimate(x0, x1, y0);
    }

    DataType id_estimate(DataType x0, DataType x1, DataType y0)
    {
      return y0;
    }
    
    DataType linear_estimate(DataType x0, DataType x1, DataType y0)
    {
      if(y0 == 0)
        return 0;
      auto sinh = is * (oldexpy1 - oldinvexpy1);
      return (y0 + A * (2 * x1 - y0) - B * sinh) / (1 + A + B * sinh / y0);
    }
    
    DataType affine_estimate(DataType x0, DataType x1, DataType y0)
    {
      auto sinh = is * (expdiode_y1_p - expdiode_y1_m);
      auto cosh = is * (expdiode_y1_p + expdiode_y1_m);

      return (x1 + R * ieq - R * sinh + R * y0 * cosh / vt) / (1 + R * C + R * cosh / vt);
    }
  };
  
  template <typename DataType>
  DiodeClipperFilterv6<DataType>::DiodeClipperFilterv6()
  :TypedBaseFilter<DataType>(1, 1)
  {
    input_delay = 1;
    output_delay = 1;
  }

  template <typename DataType>
  DiodeClipperFilterv6<DataType>::~DiodeClipperFilterv6()
  {
  }

  template <typename DataType>
  void DiodeClipperFilterv6<DataType>::setup()
  {
    Parent::setup();
    optimizer.reset(new ATK::ScalarNewtonRaphson<SimpleOverdriveFunction>(SimpleOverdriveFunction(static_cast<DataType>(1. / input_sampling_rate),
      10000, static_cast<DataType>(22e-9), static_cast<DataType>(1e-12), static_cast<DataType>(26e-3))));
  }

  template <typename DataType>
  void DiodeClipperFilterv6<DataType>::process_impl(std::size_t size) const
  {
    const DataType* ATK_RESTRICT input = converted_inputs[0];
    DataType* ATK_RESTRICT output = outputs[0];
    for(std::size_t i = 0; i < size; ++i)
    {
      optimizer->optimize(input + i, output + i);
      optimizer->get_function().update_state(input + i, output + i);
    }
  }
}

BOOST_AUTO_TEST_CASE(DiodeClipperFilterv6_const_sin1k)
{
	ATK::SimpleSinusGeneratorFilter<double> generator;
	generator.set_frequency(100);
	generator.set_output_sampling_rate(48000);
	generator.set_output_delay(6);

	ATK::OversamplingFilter<double, ATK::Oversampling6points5order_4<double> > oversampling_filter;
	oversampling_filter.set_input_sampling_rate(48000);
	oversampling_filter.set_output_sampling_rate(48000 * 4);
	oversampling_filter.set_input_port(0, &generator, 0);
	oversampling_filter.set_output_delay(1);

	ATKADC::DiodeClipperFilterv6<double> filter;
	filter.set_input_sampling_rate(48000 * 4);
	filter.set_input_port(0, &oversampling_filter, 0);
	filter.set_output_delay(5);

	ATK::IIRFilter<ATK::ButterworthLowPassCoefficients<double> > lowfilter;
	lowfilter.set_input_sampling_rate(48000 * 4);
	lowfilter.set_output_sampling_rate(48000 * 4);
	lowfilter.set_cut_frequency(48000);
	lowfilter.set_order(5);
	lowfilter.set_input_port(0, &filter, 0);

	ATK::DecimationFilter<double> decimation_filter;
	decimation_filter.set_input_sampling_rate(48000 * 4);
	decimation_filter.set_output_sampling_rate(48000);
	decimation_filter.set_input_port(0, &lowfilter, 0);

	for(size_t i = 0; i < PROCESSSIZE; ++i)
		decimation_filter.process(PROCESSSIZE);
}
