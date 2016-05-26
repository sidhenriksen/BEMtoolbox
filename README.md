This provides the development version of BEMtoolbox. There is no stable version.
BEMunit provides the base class for simulating binocular energy model units. This class also allows for extending the standard BEM framework to incorporate other nonlinearities. The demos directory gives some demonstrations for how
to use the toolbox, including how to compute a disparity tuning curve.

The documentation for the toolbox is an on-going project. If there is anything
you need that is not adequately documented or you would like to see implemented, please do get in touch [sid(dot)henriksen(at)gmail(dot)com]. 

# The components of BEMtoolbox
There are two key classes in BEMtoolbox: BEMunits and StimulusGenerators.

## BEMunit
BEMunits are your model objects. An object is created by calling:
```
bem = BEMunit();
```

You can set a bunch of model properties, such as RF size, vertical and
horizontal disparity, temporal kernel width (and type), and so on.

## StimulusGenerator
A StimulusGenerator is a superclass that has certain methods associated
with it. Crucially, each StimulusGenerator subclass *needs* a method called
generate which will return an example stereogram. For example:
```
rds = RDS()
[L,R] = rds.generate()
```

L and R is the left and right eye's image, respectively. 

It is strongly recommended (though not strictly necessary) that you inherit 
from StimulusGenerator when you create your stimulus class (see RDS.m for an 
example). The advantage of this is that there are additional functionalities 
in the StimulusGenerator class (such as easy displaying of stimuli, unit 
testing, and so on). In principle, a StimulusGenerator can be anything:
sine waves with random phases or frequencies, natural images, or, as used
here, random dot stereograms.

## Simulating responses
We do this by putting together BEMunit and StimulusGenerator
```
bem = BEMunit('x0',0,'y0',0) # sets RF centres in the middle of the stimulus
rds = RDS()
N = 1000
C = bem.simulate_spatial(rds,N) # compute responses to N stimuli (returned by calling generate on the RDS object)
```
See the demos directory for an example disparity tuning curve.

BEMtoolbox also supports spatiotemporal responses. Spatiotemporal simulations
can be done fairly easily:
```
n_frames = 20;
duration = 0.5;
bem.temporal_kernel = 'gaussian';
C = bem.simulate_spatiotemporal(rds,n_frames,duration);
```
C will here be an array which contains the temporal response of the RF to a sequence of 
generated stimuli.

Currently, there is only support for separable spacetime receptive fields
with more than two dimensions (although you can in principle use the 2D spatial RF as a 
spacetime RF with one space dimension; I might get around to implementing more proper 
support for spacetime inseparable RFs).

# Custom models
BEMtoolbox provides a flexible and extensible framework for simulating binocular
cells.

## Adding subunits
In order to add subunits, simply call the add_subunit method on the BEMunit object. E.g.
```
 bem = BEMunit();
 bem = bem.add_subunit('L_phi',pi/2,'R_phi',3/4*pi')
```
will add a subunit with a Gabor with pi/2 phase in the left and 3/4*pi phase 
in the right eye.

## Modifying subunits
All properties of the subunits are found in bem.subunits.
For example:
`bem.subunits(1).rf_params` contains the receptive field parameters for each subunit. Monocular parameters can be accessed in `bem.subunits(1).rf_params.left` and `bem.subunits(1).rf_params.right`.

One way to change the frequency of the monocular Gabors is to set these fields
and update the BEMunit object. However, this can be awkward, so a short hand
is to use:
```
subunit = 1;
bem = bem.modify_subunits(subunit,'L_f',1.6,'R_f',1.6);
```
`'L_f`' tells BEMtoolbox that we want to set the frequency (`f`) field in the 
left eye, and likewise for `'R_f'`. Frequency units are in cycles per degree.
The variable subunit indexes the subunit we want to change (in this case it is the first subunit). 

## Changing nonlinearities
BEMtoolbox supports nonlinearities at two stages of the model hierarchy:

### Simple nonlinearities 
These are the NLs that act on the simple cells in the classical
energy model framework. This is simply a squaring in the BEM, S1 = (L1+R1)^2,
where L1 and R1 refer to the responses of the monocular subunits. If we want
to change this, we simply set the `bem.subunits(i).NL` field to whatever
nonlinearity we want using a function handle, e.g.:
```
myNL = @(x)((x>0).*x); #threshold nonlinearity
bem.subunits(1).NL = myNL;
bem.subunits(2).NL = myNL;
```

### Complex nonlinearities
These are NLs that act on the complex cell response. In the original BEM, 
this is simply the identity function f(x)=x. However, multiple papers have 
explored alternatives to this (Read et al., 2002; Doi & Fujita, 2014; 
Henriksen et al., 2016, and many more). For example:
```
myNL = @(x)(x.^2);
bem.outputNL = myNL;
```
Will square the output of the energy model units (explored in 
Read et al., 2002; Henriksen et al., 2016). Alternatively,
```
theta=0.5;
gamma = 2;
myNL = @(x)( (x * (x > theta)).^gamma)
bem.outputNL = myNL
```
Will threshold and raise to power gamma.

# Advanced (slightly experimental) capabilities
## Bootstrap mode
Bootstrap mode is a way of computing responses through resampling. This is, to my knowledge,
the fastest way of running spatiotemporal simulations with spacetime separable receptive fields.
Bootstrap mode works by pre-computing the spatial responses of cells and storing these to disk.
When you then run a spatiotemporal simulation with a particular stimulus, BEMtoolbox will read the
data from disk to generate your random sequence of RDSs. It is called bootstrap mode because what
we are essentially doing is creating a bootstrap distribution of monocular spatial responses,
and resampling from this distribution before doing our temporal filtering. This is *much* faster than
generating a new RDS, computing the spatial response of each monocular subunit, and then doing the
temporal convolutions. The temporal convolutions are done through FFTs, which are very fast (nlogn). Thus, the
bottleneck is in generating the RDS and doing spatial filtering. Bootstrap mode precomputes the monocular
responses so that you only have to go through the bottleneck once per stimulus condition. Note, it is best
to use a large number of frames for generating your bootstrap distribution (ideally > 5000, depending on 
your purposes)

How it works:
```
bootstrap_mode = 1;
N = 1e3; % number of RDSs to pre-compute (large numbers are usually better)
bem = BEMunit();
rds = RDS();
C=bem.simulate_spatial(rds,N,bootstrap_mode);

```

This will save the monocular convolution values for this stimulus to disk. You will need to load the values from
disk with `load_bootstrap`; loading this explicitly means that you only have to read the values from disk once. 

If you now want to simply compute the model responses using the saved values you can do:
```
bem = bem.load_bootstrap(rds); # load bootstrapped responses
N = 0;
bootstrap_mode = 1;
C = bem.simulate_spatial(rds,N,bootstrap_mode);
```
Setting `N=0` and `bootstrap_mode=1` tells BEMtoolbox to load the responses, rather than computing new ones.
Using N>0 will make BEMtoolbox compute new monocular convolutions and *append* them to existing ones.

You can also do a spatiotemporal simulation with and without bootstrap mode. This is where bootstrap mode
is perhaps the most useful:
```
bem.temporal_kernel = 'gamma-cosine';
% 10 frames in 500 ms is 20Hz RDS.
n_frames = 10;
duration = 0.5;

# without bootstrap mode
tic
bootstrap_mode = 0; 
for i = 1:100;
    C = bem.simulate_spatiotemporal(rds,n_frames,duration,bootstrap_mode);
end
toc

# with bootstrap mode
tic
bootstrap_mode = 1; 
bem = bem.load_bootstrap(rds); # this loads the data so that we only read from disk once

for i = 1:100;
   C = bem.simulate_spatiotemporal(rds,n_frames,duration,bootstrap_mode);
end
toc
``` 
This improves the speed by more than a factor of 10.

## Parallel mode
Parallel mode allows for running `simulate_spatial` in parallel. Because
`simulate_spatiotemporal` returns an array of temporal response, there is no explicit
run_parallel parameter for this method. However, you can run `simulate_spatiotemporal` in a 
parfor loop. To run `simulate_spatial` in parallel:
```
run_parallel = 0;
C = bem.simulate_spatial(rds,N,bootstrap_mode,run_parallel);
```
As you can see, run_parallel can be used in conjunction with bootstrap_mode. This will generate
monocular responses in parallel, and then save these to disk for easy reuse (e.g. for
temporal convolutions).
