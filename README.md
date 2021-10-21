# Transfer-entropy-spectrum-in-the-Fourier-domain
This is the algorithm for calculating the transfer entropy spectrum in the Fourier-domain, which is a novel generalization of transfer entropy. There are **4 folders** in our open-source codes.

**On 2021.10.11**, we will let the MATLAB implementation of our algorithm be open source. One can learn about this implementation according to our paper.
You can find this version in the **Initial-MATLAB-version folder**. This is the algorithm implemented in our paper entitled **"Fourier-domain transfer entropy spectrum"**. You can find the pre-print of this paper via **[arXiv of Fourier-domain transfer entropy spectrum](https://arxiv.org/abs/2110.06480)**.

**On 2021.10.21**, we will relase an engineering version of the MATLAB implementation. In this version, we will accelerate the computation by matrix organization and approximation in our codes. This version is better suited for engineering using (e.g., processing videos or multiple time series) rather than scientific using. Please note that  videos are expected to be low-resolution becasue the number of time series will be a main factor of memory consumption (you can run our code no matter how many time series there are, but the large size of output reuslts will easily crash your computer). It is not completely equivalent to the approach described in the paper **"Fourier-domain transfer entropy spectrum"**. You can find this version in the **Engineering-MATLAB-version folder**. 

**On 2021.10.28**, we will release the Python version of our algorithm. You can find this version in the **Initial-Pyton-version folder**. This version is equivalent to the
the version in the **Initial-MATLAB-version branch** and designed for Python users. This is the algorithm implemented in our paper entitled **"Fourier-domain transfer entropy spectrum"**. You can find the pre-print of this paper via **[arXiv of Fourier-domain transfer entropy spectrum](https://arxiv.org/abs/2110.06480)**.

**On 2021.10.28**, we will relase an engineering version of the Python implementation. This version is equivalent to the **Engineering-MATLAB-version folder** and better suited for engineering using rather than scientific using. It is not completely equivalent to the approach described in the paper **"Fourier-domain transfer entropy spectrum"**. You can find this version in the **Engineering-Python-version folder**. 
______________________________________________________________________________________________________________________________________________________________________________

Below, we introduce the details of the **Initial-MATLAB-version folder**. Among these files, "MainFunctionforFDTES.m" 
is the main function to run. "SymbolizationFunction.m", "SearchHistory.m", "FourierDomainTransferEntropySpectrum.m", "SignificanceTestFunction.m" 
are ancillary functions to realize different steps of our approach. 

Please note that there is a file **"permutationTest.m"** uploaded as a ancillary function of "SignificanceTestFunction.m". It is not our original work, 
and we upload it only for the convenience of other users. **The credits of this code should belong to Laurens R Krol** via https://github.com/lrkrol/permutationTest. 
**When you use our algorithm, please cite Laurens R Krol's work as well**.

The function is given in following forms:
```
[TransferEntropyM,XTSymbolMatrix,XFSymbolMatrix,YTSymbolMatrix,YFSymbolMatrix,f,varargout]=MainFunctionforFDTES(X,Y,TOrder,FOrder,TLag,Type,varargin)
```
To use our code (both MATLAB and Python version), one may consider following tips:

<big>**Input:**</big>

**(1)** X is a n*1 time series vector. It can not contain Inf and NaN so you
have to pre-process it;

**(2)** Y is a n*1 time series vector. It can not contain Inf and NaN so you
have to pre-process it;

**(3)** TOrder is \lambda in the paper. It determines the coarse-graining on
time-line during the symbolization; It should be a positive number and 
never be larger than the length of time series.

**(4)** FOrder is \theta in the paper. It determines the coarse-graining on
frequency-line during the symbolization; It should be a positive number  
and never be larger than the length of frequency bands.

**(5)** TLag is \beta in the paper. It determines the length historical 
information added in transfer entropy calculation;

**(6)** Type is used to determine if you want to set all negative TEs in the
TransferEntropyM as 0. If Type=1, then all negative TEs become 0. If
Type=0, then we keep the raw data and you can still obtain negative TEs.

**(7)** varargin contains the information of statistic significance test. If
you do not want the statistic significance test, just do not enter
varargin. If varargin is not empty, it is expected to contain two inputs
so that NumofSurrogates=varargin{1} and ShiftInterval=varargin{2}. NumofSurrogates 
determines how many surrogates you are going to generate during the 
statistic significance test. Please note that surrogate generation and 
statistic significance test is computationally costly. For instance, 
NumofSurrogates=k means that k surrogates are generated, k times of 
Fourier-domain transfer entropy spectrum calculations are implemented, 
and a permutation test is carried out. ShiftInterval should be a vector 
[p,q], where p is the minimum time shift and q is the maximum time shift. 

<big>**Output:**</big>

**(1)** TransferEntropyM is the Fourier-domain transfer entropy spectrum
T(X,Y,t,\omega). If you calculate sum(TransferEntropyM), you obtain
T(X,Y,t); If you calculate sum(TransferEntropyM,2), you obtain
T(X,Y,\omega); If you calculate sum(sum(TransferEntropyM)), you obtain
T(X,Y); 

**(2)** XTSymbolMatrix is the symbolization result of X along the time line,
namely \Lambda_{\mathsf{X}}\left(t,\omega\right) in the paper.

**(3)** XFSymbolMatrix is the symbolization result of X along the frequency 
line, namely \Theta_{\mathsf{X}}\left(t,\omega\right) in the paper.

**(4)** YTSymbolMatrix is the symbolization result of Y along the time line,
namely \Theta_{\mathsf{Y}}\left(t,\omega\right) in the paper.

**(5)** YFSymbolMatrix is the symbolization result of Y along the frequency 
line, namely \Theta_{\mathsf{Y}}\left(t,\omega\right) in the paper.

**(6)** f is the frequency band obtained from wavelet transfrom. Please note
that you are not required to specify a sampling frequency so f is
returned as number of cycles per sample. f is \omega in the paper.

**(7)** If you want to do statistic significance test and you do provide the
varargin, then you can ask for varargout. varargout{1}=P denotes the p
value of significance test. varargout{2}=Effectsize denotes the effect 
size. varargout{3}=TESample denotes the set of T(X,Y,\omega). varargout{4}
denotes the the set of T(X',Y,\omega) generated by surrogates.

<big>**Examples:**</big>
```
(1) If you want all: 
[TransferEntropyM,XTSymbolMatrix,XFSymbolMatrix,YTSymbolMatrix,YFSymbolMatrix,f,P,Effectsize,TESample,SurrogatesTESample]=MainFunctionforFDTES(X,Y,TOrder,FOrder,TLag,Type,NumofSurrogates,ShiftInterval)

(2) If you do not want the statistic significance test: 
[TransferEntropyM,XTSymbolMatrix,XFSymbolMatrix,YTSymbolMatrix,YFSymbolMatrix,f]=MainFunctionforFDTES(X,Y,TOrder,FOrder,TLag,Type)
```
______________________________________________________________________________________________________________________________________________________________________________

Below, we introduce the details of the **Engineering-MATLAB-version folder**. Among these files, "EngineeringMFforFDTES.m" 
is the main function to run. Other files are ancillary functions to realize different steps of our approach. Please note that 
the default setting of this code includes no statistic significance test to control memory consumption. Of course, you can 
add statistic significance test manually.

The function is given in following forms:
```
[ATransferEntropyM,TFCoordM,TFCMatrix]=EngineeringMFforFDTES(X,Wavelet,DLevel,TOrder,FOrder,SymbolizationType,SamplingSize,TLag,Type)
```

To use our code (both MATLAB and Python version), one may consider following tips:

<big>**Input:**</big>

**(1)** X is a h*k*n time series matrix. It can not contain Inf and NaN so you
have to pre-process it. n is the length of time series and h*k is the number
of time series. In practice, (a) if you want to process a video or multiple time series, 
then you can try this code and j*k is the number of pixels or single
time series; (b) however, you are not recommend to use this code if X
contains too many time series, e.g., j*k>1000. This code can still run
in such situations but it will cost unacceptable computer memory. In such a case, you 
can try the simplified algorithm version;

**(2)** Wavelet is the name of wavelet function. For instance, you can input
Wavelet='sym4' for a near symmetric wavelet. In the document of function
wfilters, you can see the built in wavelet functions in MATLAB. Please
use one of them according to your demands.

**(3)** DLevel is the number of levels in wavelet decomposition.

**(4)** TOrder is \lambda in the paper. It determines the coarse-graining on
time-line during the symbolization; It should be a positive number and
never be larger than the length of time series.

**(5)** FOrder is \theta in the paper. It determines the coarse-graining on
frequency-line during the symbolization; It should be a positive number
and never be larger than the length of frequency bands. In this version,
you must let FOrder<DLevel+1. In practice, we suggest that you can
choose FOrder<0.5*(DLevel+1)

**(6)** SymbolizationType is a variable to determine the type of symbolization
approach. You can select SymbolizationType=1 for the standard symbolization
or SymbolizationType=1 for the simplified symbolization

**(7)** SamplingSize is the sampling size used in random calculation of 
Fourier-domain transfer entropy spectrum

**(8)** TLag is \beta in the paper. It determines the length historical
information added in transfer entropy calculation;

**(9)** Type is used to determine if you want to set all negative TEs in the
TransferEntropyM as 0. If Type=1, then all negative TEs become 0. If
Type=0, then we keep the raw data and you can still obtain negative TEs.

<big>**Output:**</big>

**(1)** ATransferEntropyM is the (h*k)*(h*k)*a*b matrix (it is four dimensional!!)
of the Fourier-domain transfer entropy spectrum T(X(i_{1},j_{1},:),X(i_{2},j_{2},:),t,\omega).
If you search across the first and two dimension, you can find the Fourier-domain transfer entropy
spectrum T(X(i_{1},j_{1},:),X(i_{2},j_{2},:),t,\omega) between any two time series X(i_{1},j_{1},:)
and X(i_{2},j_{2},:) at the location "location=sub2ind([(h*k),(h*k)],i,j)". Therefore, you
can get the spectrum by TransferEntropyM(:,:,sub2ind([(h*k),(h*k)],i,j)). Here a
and b denotes the length of frequency and time line, respectively.

**(2)** TFCoordM is a (h*k)*(h*k)*a*b matrix (it is four dimensional!!) of the
symbolization result along the time and frequency lines. XTSymbolMatrix(i,j,:,:)
is the symbolization of a variable (i,j).

**(3)** TFCMatrix is a (h*k)*(h*k)*a*b matrix of the linear correlation result of X.
