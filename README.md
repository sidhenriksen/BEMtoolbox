This provides the development version of BEMtoolbox. BEMunit provides the
base class for simulating binocular energy model units. This class also 
allows for extending the standard BEM framework to incorporate arbitrary
nonlinearities. The demos directory gives some demonstrations for how
to use the toolbox, including how to compute a disparity tuning curve.

# The components of BEMtoolbox
There are two key classes in BEMtoolbox: BEMunits and StimulusGenerators.

## BEMunit
BEMunits are your model objects. An object is created by calling:
```
bem = BEMunit();


You can set a bunch of model properties, such as RF size, vertical and
horizontal disparity, temporal kernel width (and type), and so on.

## StimulusGenerator
A StimulusGenerator is a superclass that has certain methods associated
with it. Crucially, each StimulusGenerator subclass *needs* a method called
generate which will return an example stereogram. For example:
```
rds = RDS()
[L,R] = rds.generate()

L and R is the left and right eye's image, respectively. 

It is strongly recommended (though not strictly necessary) that you inherit 
from StimulusGenerator when you create your stimulus class (see RDS.m for an 
example). The advantage of this is that there are additional functionalities 
in the StimulusGenerator class (such as easy displaying of stimuli, unit 
testing, and so on).

## Putting together BEMunit and StimulusGenerator
```
 bem = BEMunit('x0',0,'y0',0) # sets RF centres in the middle of the stimulus
 rds = RDS()
 N = 1000
 C = bem.simulate_spatial(rds,N) # compute responses to N stimuli (returned by calling generate on the RDS object)

See the demos directory for an example disparity tuning curve.

# Custom models
BEMtoolbox provides a flexible and extensible framework for simulating binocular
cells.

## Adding subunits
In order to add subunits, simply call the add_subunit method on the BEMunit object. E.g.
```
 bem = BEMunit();
 bem = bem.add_subunit('L_phi',pi/2,'R_phi',3/4*pi')
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

'L_f' tells BEMtoolbox that we want to set the frequency (f) field in the 
left eye, and likewise for 'R_f'. The variable subunit indexes the subunit
we want to change (in this case it is the first subunit). 

## Changing nonlinearities
BEMtoolbox supports nonlinearities at two stages of the model hierarchy:

### Simple nonlinearities 
These are the NLs that act on the simple cells in the classical
energy model framework. This is simply a squaring in the BEM, S1 = (L1+R1)^2,
where L1 and R1 refer to the responses of the monocular subunits. If we want
to change this, we simply set the bem.subunits(i).NL field to whatever
nonlinearity we want, e.g.:
```
myNL = @(x)((x>0).*x); #threshold nonlinearity
bem.subunits(1).NL = myNL;
bem.subunits(2).NL = myNL;

### Complex nonlinearities
These are NLs that act on the complex cell response. In the original BEM, 
this is simply the identity function. However, multiple papers have explored
alternatives to this (Read et al., 2002; Doi & Fujita, 2014; 
Henriksen et al., 2016, and many more). For example:
```
myNL = @(x)(x.^2);
bem.outputNL = myNL;

Will square the output of the energy model units. Alternatively,
```
theta=0.5;
gamma = 2;
myNL = @(x)( (x * (x > theta)).^gamma)
bem.outputNL = myNL

Will threshold and raise to power gamma.