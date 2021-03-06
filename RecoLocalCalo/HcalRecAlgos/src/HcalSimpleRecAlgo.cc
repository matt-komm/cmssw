#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSimpleRecAlgo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"
#include <algorithm> // for "max"
#include <math.h>

static double MaximumFractionalError = 0.0005; // 0.05% error allowed from this source

HcalSimpleRecAlgo::HcalSimpleRecAlgo(int firstSample, int samplesToAdd, bool correctForTimeslew, bool correctForPulse, float phaseNS) : 
  firstSample_(firstSample), 
  samplesToAdd_(samplesToAdd), 
  correctForTimeslew_(correctForTimeslew) {
  if (correctForPulse) 
    pulseCorr_=std::auto_ptr<HcalPulseContainmentCorrection>(new HcalPulseContainmentCorrection(samplesToAdd_,phaseNS,MaximumFractionalError));
}

HcalSimpleRecAlgo::HcalSimpleRecAlgo(int firstSample, int samplesToAdd) : 
  firstSample_(firstSample), 
  samplesToAdd_(samplesToAdd), 
  correctForTimeslew_(false) {
}


void HcalSimpleRecAlgo::resetTimeSamples(int firstSample, int samplesToAdd)
{
  firstSample_ = firstSample;
  samplesToAdd_ = samplesToAdd;
}


///Timeshift correction for HPDs based on the position of the peak ADC measurement.
///  Allows for an accurate determination of the relative phase of the pulse shape from
///  the HPD.  Calculated based on a weighted sum of the -1,0,+1 samples relative to the peak
///  as follows:  wpksamp = (0*sample[0] + 1*sample[1] + 2*sample[2]) / (sample[0] + sample[1] + sample[2])
///  where sample[1] is the maximum ADC sample value.
static float timeshift_ns_hbheho(float wpksamp);

///Same as above, but for the HF PMTs.
static float timeshift_ns_hf(float wpksamp);


namespace HcalSimpleRecAlgoImpl {
  template<class Digi, class RecHit>
  inline RecHit reco(const Digi& digi, const HcalCoder& coder, const HcalCalibrations& calibs, 
		     int ifirst, int n, bool slewCorrect, const HcalPulseContainmentCorrection* corr,
		     HcalTimeSlew::BiasSetting slewFlavor) {
    CaloSamples tool;
    coder.adc2fC(digi,tool);

    double ampl=0; int maxI = -1; double maxA = -1e10; float ta=0;
    double fc_ampl=0;
    for (int i=ifirst; i<tool.size() && i<n+ifirst; i++) {
      int capid=digi[i].capid();
      ta = (tool[i]-calibs.pedestal(capid)); // pedestal subtraction
      fc_ampl+=ta; 
      ta*= calibs.respcorrgain(capid) ; // fC --> GeV
      ampl+=ta;
      if(ta>maxA){
	maxA=ta;
	maxI=i;
      }
    }

    float time = -9999;
    ////Cannot calculate time value with max ADC sample at first or last position in window....
    if(maxI==0 || maxI==(tool.size()-1)) {      
      LogDebug("HCAL Pulse") << "HcalSimpleRecAlgo::reconstruct :" 
					       << " Invalid max amplitude position, " 
					       << " max Amplitude: "<< maxI
					       << " first: "<<ifirst
					       << " last: "<<(tool.size()-1)
					       << std::endl;
    } else {
      int capid=digi[maxI-1].capid();
      float t0 = ((tool[maxI-1]-calibs.pedestal(capid))*calibs.respcorrgain(capid) );
      capid=digi[maxI+1].capid();
      float t2 = ((tool[maxI+1]-calibs.pedestal(capid))*calibs.respcorrgain(capid) );

      // Handle negative excursions by moving "zero":
      float minA=t0;
      if (maxA<minA) minA=maxA;
      if (t2<minA)   minA=t2;
      if (minA<0) { maxA-=minA; t0-=minA; t2-=minA; } // positivizes all samples

      float wpksamp = (t0 + maxA + t2);
      if (wpksamp!=0) wpksamp=(maxA + 2.0*t2) / wpksamp; 
      time = (maxI - digi.presamples())*25.0 + timeshift_ns_hbheho(wpksamp);
      if (corr!=0) {
	// Apply phase-based amplitude correction:
	ampl *= corr->getCorrection(fc_ampl);
	//      std::cout << fc_ampl << " --> " << corr->getCorrection(fc_ampl) << std::endl;
      }
      if (slewCorrect) time-=HcalTimeSlew::delay(std::max(1.0,fc_ampl),slewFlavor);

      time=time-calibs.timecorr(); // time calibration
    }

    return RecHit(digi.id(),ampl,time);    
  }
}

HBHERecHit HcalSimpleRecAlgo::reconstruct(const HBHEDataFrame& digi, const HcalCoder& coder, const HcalCalibrations& calibs) const {
  return HcalSimpleRecAlgoImpl::reco<HBHEDataFrame,HBHERecHit>(digi,coder,calibs,
							       firstSample_,samplesToAdd_,correctForTimeslew_,
							       pulseCorr_.get(),
							       HcalTimeSlew::Medium);
}

HORecHit HcalSimpleRecAlgo::reconstruct(const HODataFrame& digi, const HcalCoder& coder, const HcalCalibrations& calibs) const {
  return HcalSimpleRecAlgoImpl::reco<HODataFrame,HORecHit>(digi,coder,calibs,
							   firstSample_,samplesToAdd_,correctForTimeslew_,
							   pulseCorr_.get(),
							   HcalTimeSlew::Slow);
}

HcalCalibRecHit HcalSimpleRecAlgo::reconstruct(const HcalCalibDataFrame& digi, const HcalCoder& coder, const HcalCalibrations& calibs) const {
  return HcalSimpleRecAlgoImpl::reco<HcalCalibDataFrame,HcalCalibRecHit>(digi,coder,calibs,
									 firstSample_,samplesToAdd_,correctForTimeslew_,
									 pulseCorr_.get(),
									 HcalTimeSlew::Fast);
}

HFRecHit HcalSimpleRecAlgo::reconstruct(const HFDataFrame& digi, const HcalCoder& coder, const HcalCalibrations& calibs) const {
  CaloSamples tool;
  coder.adc2fC(digi,tool);

  double ampl=0; int maxI = -1; double maxA = -1e10; float ta=0; float amp_fC=0;
  for (int i=firstSample_; i<tool.size() && i<samplesToAdd_+firstSample_; i++) {
    int capid=digi[i].capid();
    ta = (tool[i]-calibs.pedestal(capid))*calibs.respcorrgain(capid);
    ampl+=ta;
    amp_fC += tool[i]-calibs.pedestal(capid);
    if(ta>maxA){
      maxA=ta;
      maxI=i;
    }
  }

  float time=-9999.0;
  ////Cannot calculate time value with max ADC sample at first or last position in window....
  if(maxI==0 || maxI==(tool.size()-1)) {
      LogDebug("HCAL Pulse") << "HcalSimpleRecAlgo::reconstruct :" 
					       << " Invalid max amplitude position, " 
					       << " max Amplitude: "<< maxI
					       << " first: "<<firstSample_
					       << " last: "<<(tool.size()-1)
					       << std::endl;
  } else {
    int capid=digi[maxI-1].capid();
    float t0 = (tool[maxI-1]-calibs.pedestal(capid))*calibs.respcorrgain(capid);
    capid=digi[maxI+1].capid();
    float t2 = (tool[maxI+1]-calibs.pedestal(capid))*calibs.respcorrgain(capid);

    // Handle negative excursions by moving "zero":
    float zerocorr=std::min(t0,t2);
    if (zerocorr<0.f) {
      t0   -= zerocorr;
      t2   -= zerocorr;
      maxA -= zerocorr;
    }
    
    // pair the peak with the larger of the two neighboring time samples
    float wpksamp=0.f;
    if (t0>t2) {
      wpksamp = t0+maxA;
      if (wpksamp != 0.f) wpksamp = maxA/wpksamp;
    } else {
      wpksamp = maxA+t2;
      if (wpksamp != 0.f) wpksamp = 1.+(t2/wpksamp);
    }

    time = (maxI - digi.presamples())*25.0 + timeshift_ns_hf(wpksamp);

    if (correctForTimeslew_ && (amp_fC>0)) {
      // -5.12327 - put in calibs.timecorr()
      double tslew=exp(0.337681-5.94689e-4*amp_fC)+exp(2.44628-1.34888e-2*amp_fC);
      time -= (float)tslew;
    }

    time=time-calibs.timecorr();
  }

  return HFRecHit(digi.id(),ampl,time); 
}

// timeshift implementation

static const float wpksamp0_hbheho = 0.5;
static const int   num_bins_hbheho = 61;

static const float actual_ns_hbheho[num_bins_hbheho] = {
-5.44000, // 0.500, 0.000-0.017
-4.84250, // 0.517, 0.017-0.033
-4.26500, // 0.533, 0.033-0.050
-3.71000, // 0.550, 0.050-0.067
-3.18000, // 0.567, 0.067-0.083
-2.66250, // 0.583, 0.083-0.100
-2.17250, // 0.600, 0.100-0.117
-1.69000, // 0.617, 0.117-0.133
-1.23000, // 0.633, 0.133-0.150
-0.78000, // 0.650, 0.150-0.167
-0.34250, // 0.667, 0.167-0.183
 0.08250, // 0.683, 0.183-0.200
 0.50250, // 0.700, 0.200-0.217
 0.90500, // 0.717, 0.217-0.233
 1.30500, // 0.733, 0.233-0.250
 1.69500, // 0.750, 0.250-0.267
 2.07750, // 0.767, 0.267-0.283
 2.45750, // 0.783, 0.283-0.300
 2.82500, // 0.800, 0.300-0.317
 3.19250, // 0.817, 0.317-0.333
 3.55750, // 0.833, 0.333-0.350
 3.91750, // 0.850, 0.350-0.367
 4.27500, // 0.867, 0.367-0.383
 4.63000, // 0.883, 0.383-0.400
 4.98500, // 0.900, 0.400-0.417
 5.33750, // 0.917, 0.417-0.433
 5.69500, // 0.933, 0.433-0.450
 6.05000, // 0.950, 0.450-0.467
 6.40500, // 0.967, 0.467-0.483
 6.77000, // 0.983, 0.483-0.500
 7.13500, // 1.000, 0.500-0.517
 7.50000, // 1.017, 0.517-0.533
 7.88250, // 1.033, 0.533-0.550
 8.26500, // 1.050, 0.550-0.567
 8.66000, // 1.067, 0.567-0.583
 9.07000, // 1.083, 0.583-0.600
 9.48250, // 1.100, 0.600-0.617
 9.92750, // 1.117, 0.617-0.633
10.37750, // 1.133, 0.633-0.650
10.87500, // 1.150, 0.650-0.667
11.38000, // 1.167, 0.667-0.683
11.95250, // 1.183, 0.683-0.700
12.55000, // 1.200, 0.700-0.717
13.22750, // 1.217, 0.717-0.733
13.98500, // 1.233, 0.733-0.750
14.81500, // 1.250, 0.750-0.767
15.71500, // 1.267, 0.767-0.783
16.63750, // 1.283, 0.783-0.800
17.53750, // 1.300, 0.800-0.817
18.38500, // 1.317, 0.817-0.833
19.16500, // 1.333, 0.833-0.850
19.89750, // 1.350, 0.850-0.867
20.59250, // 1.367, 0.867-0.883
21.24250, // 1.383, 0.883-0.900
21.85250, // 1.400, 0.900-0.917
22.44500, // 1.417, 0.917-0.933
22.99500, // 1.433, 0.933-0.950
23.53250, // 1.450, 0.950-0.967
24.03750, // 1.467, 0.967-0.983
24.53250, // 1.483, 0.983-1.000
25.00000  // 1.500, 1.000-1.017 - keep for interpolation
};

float timeshift_ns_hbheho(float wpksamp) {
  float flx = (num_bins_hbheho-1)*(wpksamp - wpksamp0_hbheho);
  int index = (int)flx;
  float yval;

  if      (index <    0)               return actual_ns_hbheho[0];
  else if (index >= num_bins_hbheho-1) return actual_ns_hbheho[num_bins_hbheho-1];

  // else interpolate:
  float y1 = actual_ns_hbheho[index];
  float y2 = actual_ns_hbheho[index+1];

  yval = y1 + (y2-y1)*(flx-(float)index);

  return yval;
}

static const int   num_bins_hf = 101;
static const float wpksamp0_hf = 0.5;

static const float actual_ns_hf[num_bins_hf] = {
 0.00250, // 0.000-0.010
 0.04500, // 0.010-0.020
 0.08750, // 0.020-0.030
 0.13000, // 0.030-0.040
 0.17250, // 0.040-0.050
 0.21500, // 0.050-0.060
 0.26000, // 0.060-0.070
 0.30250, // 0.070-0.080
 0.34500, // 0.080-0.090
 0.38750, // 0.090-0.100
 0.42750, // 0.100-0.110
 0.46000, // 0.110-0.120
 0.49250, // 0.120-0.130
 0.52500, // 0.130-0.140
 0.55750, // 0.140-0.150
 0.59000, // 0.150-0.160
 0.62250, // 0.160-0.170
 0.65500, // 0.170-0.180
 0.68750, // 0.180-0.190
 0.72000, // 0.190-0.200
 0.75250, // 0.200-0.210
 0.78500, // 0.210-0.220
 0.81750, // 0.220-0.230
 0.85000, // 0.230-0.240
 0.88250, // 0.240-0.250
 0.91500, // 0.250-0.260
 0.95500, // 0.260-0.270
 0.99250, // 0.270-0.280
 1.03250, // 0.280-0.290
 1.07000, // 0.290-0.300
 1.10750, // 0.300-0.310
 1.14750, // 0.310-0.320
 1.18500, // 0.320-0.330
 1.22500, // 0.330-0.340
 1.26250, // 0.340-0.350
 1.30000, // 0.350-0.360
 1.34000, // 0.360-0.370
 1.37750, // 0.370-0.380
 1.41750, // 0.380-0.390
 1.48750, // 0.390-0.400
 1.55750, // 0.400-0.410
 1.62750, // 0.410-0.420
 1.69750, // 0.420-0.430
 1.76750, // 0.430-0.440
 1.83750, // 0.440-0.450
 1.90750, // 0.450-0.460
 2.06750, // 0.460-0.470
 2.23250, // 0.470-0.480
 2.40000, // 0.480-0.490
 2.82250, // 0.490-0.500
 3.81000, // 0.500-0.510
 6.90500, // 0.510-0.520
 8.99250, // 0.520-0.530
10.50000, // 0.530-0.540
11.68250, // 0.540-0.550
12.66250, // 0.550-0.560
13.50250, // 0.560-0.570
14.23750, // 0.570-0.580
14.89750, // 0.580-0.590
15.49000, // 0.590-0.600
16.03250, // 0.600-0.610
16.53250, // 0.610-0.620
17.00000, // 0.620-0.630
17.44000, // 0.630-0.640
17.85250, // 0.640-0.650
18.24000, // 0.650-0.660
18.61000, // 0.660-0.670
18.96750, // 0.670-0.680
19.30500, // 0.680-0.690
19.63000, // 0.690-0.700
19.94500, // 0.700-0.710
20.24500, // 0.710-0.720
20.54000, // 0.720-0.730
20.82250, // 0.730-0.740
21.09750, // 0.740-0.750
21.37000, // 0.750-0.760
21.62750, // 0.760-0.770
21.88500, // 0.770-0.780
22.13000, // 0.780-0.790
22.37250, // 0.790-0.800
22.60250, // 0.800-0.810
22.83000, // 0.810-0.820
23.04250, // 0.820-0.830
23.24500, // 0.830-0.840
23.44250, // 0.840-0.850
23.61000, // 0.850-0.860
23.77750, // 0.860-0.870
23.93500, // 0.870-0.880
24.05500, // 0.880-0.890
24.17250, // 0.890-0.900
24.29000, // 0.900-0.910
24.40750, // 0.910-0.920
24.48250, // 0.920-0.930
24.55500, // 0.930-0.940
24.62500, // 0.940-0.950
24.69750, // 0.950-0.960
24.77000, // 0.960-0.970
24.84000, // 0.970-0.980
24.91250, // 0.980-0.990
24.95500, // 0.990-1.000
24.99750, // 1.000-1.010 - keep for interpolation
};

float timeshift_ns_hf(float wpksamp) {
  float flx = (num_bins_hf-1)*(wpksamp-wpksamp0_hf);
  int index = (int)flx;
  float yval;
  
  if      (index <  0)             return actual_ns_hf[0];
  else if (index >= num_bins_hf-1) return actual_ns_hf[num_bins_hf-1];

  // else interpolate:
  float y1       = actual_ns_hf[index];
  float y2       = actual_ns_hf[index+1];

  // float delta_x  = 1/(float)num_bins_hf;
  // yval = y1 + (y2-y1)*(flx-(float)index)/delta_x;

  yval = y1 + (y2-y1)*(flx-(float)index);
  return yval;
}
