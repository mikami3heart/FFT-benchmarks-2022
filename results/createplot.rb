# To use this first install the GR plotting library
# https://github.com/sciapp/gr
#
# Then install Ruby (it is likely available in your
# operating systems package managers, but installing
# from source should also work)
# https://www.ruby-lang.org/en/downloads/
# 
# Finally install GR using
# >gem install gr
#
# Set the gr directory and plotting mode
# >export GKS_WSTYPE=100
# >export GRDIR="PATH_TO_GR_INSTALLATION"
#
# Run the script to produce the plots
# >ruby createplot.rb
 
require 'gr/plot'
require 'csv'

fftw3_raw = CSV.read("serial/fftw3/1d/FFT-Fugaku-serial.FFTW3.1D.csv", headers: true)
ssl2_raw = CSV.read("serial/ssl2/1d/FFT-Fugaku-serial.SSL2.1D.csv", headers: true)
# Start with 4th endtry due to missing values and imprecise timing
# For FFTW3, largest sizes timeout in setup phase
fftw3patient_size = fftw3_raw['FFT_SIZE'][3..-3].collect { |val| val.to_f }
fftw3patient_time = fftw3_raw['EXECUTE BWD FFTW_PATIENT'][3..-3].collect { |val| val.to_f*1000.0 }
fftw3estimate_size = fftw3_raw['FFT_SIZE'][3..-1].collect { |val| val.to_f }
fftw3estimate_time = fftw3_raw['EXECUTE BWD FFT_ESTIMATE'][3..-1].collect { |val| val.to_f*1000.0 }
ssl2_size = ssl2_raw['FFT_SIZE'][3..-1].collect { |val| val.to_f }
ssl2_time = ssl2_raw['BWD'][3..-1].collect { |val| val.to_f*1000.0 }
# Start with 4th entry due to missing values and imprecise timing
GR.plot([fftw3patient_size,fftw3patient_time, spec: 'ro-'],
        [fftw3estimate_size,fftw3estimate_time, spec: 'bd-'],
        [ssl2_size,ssl2_time, spec: 'g*-'],
       title: "Comparison of FFTW3 and SSL2 1D Inverse FFT performance on Fugaku",
       xlabel: "Size of FFT (Gridpoints)",
       ylabel: "Time for FFT (ms)",
       labels: ["FFTW3 Patient","FFTW3 Estimate","SSL2"],
       xlim: [0.8*[fftw3patient_size.min,fftw3estimate_size.min,ssl2_size.min].min,
              1.2*[fftw3patient_size.max,fftw3estimate_size.max,ssl2_size.max].max],
       ylim: [0.8*[fftw3patient_time.min,fftw3estimate_time.min,ssl2_time.min].min,
              1.2*[fftw3patient_time.max,fftw3estimate_time.min,ssl2_time.max].max],
       xlog: true,
       ylog: true,
       location: 2,
)
GR.savefig("Comparison1D.svg")
