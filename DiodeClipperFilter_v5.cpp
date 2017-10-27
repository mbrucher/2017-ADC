/**
 * \file DiodeClipperFilter_v5.cpp
 */

#include "DiodeClipperFilter_v5.h"

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
  class DiodeClipperFilterv5<DataType_>::SimpleOverdriveFunction
  {
  public:
    typedef DataType_ DataType;
  protected:
    DataType A;
    DataType B;
    DataType is;
    DataType vt;

    DataType oldy0;
    DataType oldexpy0;
    DataType oldinvexpy0;
    DataType oldy1;
    DataType oldexpy1;
    DataType oldinvexpy1;
  public:
    SimpleOverdriveFunction(DataType dt, DataType R, DataType C, DataType is, DataType vt)
    :is(is), vt(vt)
    {
      A = dt / (2 * C * R);
      B = dt / (2 * C);

      oldy0 = oldy1 = 0;
      oldexpy0 = oldinvexpy0 = oldexpy1 = oldinvexpy1 = 1;
    }
    
    std::pair<DataType, DataType> operator()(const DataType* ATK_RESTRICT input, DataType* ATK_RESTRICT output, DataType y1)
    {
      auto x1 = input[0];
      auto y0 = output[-1];
      DataType expdiode_y1_p = fmath::exp(y1 / vt);
      DataType expdiode_y1_m = 1 / expdiode_y1_p;
      DataType expdiode_y0_p;
      DataType expdiode_y0_m;
	  
	    if(y0 == oldy0)
	    {
	      expdiode_y0_p = oldexpy0;
	      expdiode_y0_m = oldinvexpy0;
	    }
	    else if(y0 == oldy1)
	    {
	      expdiode_y0_p = oldexpy1;
	      expdiode_y0_m = oldinvexpy1;
	    }
	    else
	    {
	      expdiode_y0_p = fmath::exp(y0 / vt);
	      expdiode_y0_m = 1 / expdiode_y0_p;
	    }
	  
      oldy0 = y0;
      oldexpy0 = expdiode_y0_p;
      oldinvexpy0 = expdiode_y0_m;

      oldy1 = y1;
      oldexpy1 = expdiode_y1_p;
      oldinvexpy1 = expdiode_y1_m;
      
      std::pair<DataType, DataType> diode = std::make_pair(is * (expdiode_y1_p - expdiode_y1_m), is * (expdiode_y1_p + expdiode_y1_m) / vt);
      auto old_diode = is * (expdiode_y0_p - expdiode_y0_m);
      return std::make_pair(y0 - y1 + 2 * A * x1 - (A * (y1 + y0) + B * (diode.first + old_diode)), -1 - A - B * diode.second);
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
      auto sinh = is * (oldexpy1 - oldinvexpy1);
      auto cosh = is * (oldexpy1 + oldinvexpy1);
      return (y0 + A * (2 * x1 - y0) - B * (2 * sinh - y0 / vt * cosh) ) / (B * cosh / vt + 1 + A);
    }
  };
  
  template <typename DataType>
  DiodeClipperFilterv5<DataType>::DiodeClipperFilterv5()
  :TypedBaseFilter<DataType>(1, 1)
  {
    input_delay = 1;
    output_delay = 1;
  }

  template <typename DataType>
  DiodeClipperFilterv5<DataType>::~DiodeClipperFilterv5()
  {
  }

  template <typename DataType>
  void DiodeClipperFilterv5<DataType>::setup()
  {
    Parent::setup();
    optimizer.reset(new ATK::ScalarNewtonRaphson<SimpleOverdriveFunction>(SimpleOverdriveFunction(static_cast<DataType>(1. / input_sampling_rate),
      10000, static_cast<DataType>(22e-9), static_cast<DataType>(1e-12), static_cast<DataType>(26e-3))));
  }

  template <typename DataType>
  void DiodeClipperFilterv5<DataType>::process_impl(std::size_t size) const
  {
    const DataType* ATK_RESTRICT input = converted_inputs[0];
    DataType* ATK_RESTRICT output = outputs[0];
    for(std::size_t i = 0; i < size; ++i)
    {
      optimizer->optimize(input + i, output + i);
    }
  }
}

BOOST_AUTO_TEST_CASE(DiodeClipperFilterv5_const_sin1k)
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

	ATKADC::DiodeClipperFilterv5<double> filter;
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
