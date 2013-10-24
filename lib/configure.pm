package configure;

use Physics::ElectronProp::Auxiliary ':constants';

our $rad	= 1*10**-2;
our $E0  	= 1.9 * 10**5;
our $w		= 1.055*10**9;

our $cavity1 = Physics::ElectronProp::RF_Cavity->new(
  sol_name	=> 'cav',
  radius 	=> $rad,
  lens_length	=> 1*10**-2,
  front_pos	=> 0*10**-2,
  mode		=> {field => 'TM', m => 0, n => 1, p => 0},
  epsilon	=> epsilon_0,
  mu		=> mu_0,
  omega		=> $w,
  E_0		=> $E0,
);

our $cavity2 = Physics::ElectronProp::RF_Cavity->new(
  sol_name	=> 'cav',
  radius 	=> $rad,
  lens_length	=> 1*10**-2,
  front_pos	=> 3*10**-2,
  mode		=> {field => 'TM', m => 0, n => 1, p => 0},
  epsilon	=> epsilon_0,
  mu		=> mu_0,
  omega		=> $w,
  E_0		=> $E0,
);

our $electron1 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.1*10**-3,0,0*10**-2],
);

our $electron2 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.25*10**-3,0,0*10**-2],
);

our $electron3 = Physics::ElectronProp::Electron->new(
  energy 	=> 50*10**3,
  position 	=> [.5*10**-3,0,0*10**-2],
);

$cavity1->{phase} = - pi/2 + $cavity1->lens_length*$cavity1->omega/($electron1->velocity->index(2))/2;
$cavity2->{phase} = - pi/2 + $cavity2->omega*($cavity1->lens_length/2+$cavity2->front_pos)/($electron1->velocity->index(2));


our $sim = Physics::ElectronProp::Simulation_Cart->new(
  electrons	=> [$electron2,$electron1,$electron3],
  lens		=> [$cavity1,$cavity2],
  z_start	=> 0*10**-2,
  z_end		=> 7.5*10**-2,
  steps		=> 10000,
  instant	=> undef,
);
return $sim;
