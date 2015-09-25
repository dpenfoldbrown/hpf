use strict;
use Carp;
use vars qw( $AUTOLOAD $obj_table );
use DDB::UTIL;
use DDB::SEQUENCE;
use DDB::YRC;

my $id=shift;
my $width=800;
initialize_ddb();
my $seq = DDB::SEQUENCE->new(id=>$id);
$ddb_global{dbh} = DDB::UTIL->connect_db(db=>'hpf');
$seq->load();
my $sseq = $seq->get_sseq();
my $svg = DDB::YRC->_displaySequenceSvg(sseq=>$sseq, width=>$width);
print $svg
