package Physics::ElectronProp::Auxiliary;

use warnings;
use strict;

use parent 'Exporter';
our %EXPORT_TAGS = ( 
  constants   => [ qw/ pi me qe epsilon_0 vc mu_o / ],
);

## Constants ##

use constant {
  pi => 4 * atan2(1,1),
  me => 9.11e-31,
  qe => 1.6e-19,
  epsilon_0 => 8.85e-12,
  vc => 2.9979e8,
  mu_0 => 1/epsilon_0/vc**2,
};
