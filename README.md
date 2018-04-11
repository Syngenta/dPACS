# dPACS.pl


This perl script is intended as a design tool for use with the dPACS method [insert literature reference here - TBC]

### Installation

Installation assumes a basic knowledge of perl and the EMBOSS package.

In addition to this script, you will need to install and correctly configure the EMBOSS package (<http://emboss.open-bio.org/>)  along with the REBASE database (<http://rebase.neb.com/rebase/rebase.html>) which should be configured using the EMBOSS rebaseextract program.
In addition, the dPACS script should be edited (near the top) to include the location of a REBASE proto format file, which should be the same version as is imported into EMBOSS, though it can be a separate copy.

You'll need the Docopt, List::Util and IPC::Open2 perl modules installed. 

The script also requires access to temp space (e.g./tmp) if you need to change this, it's at the top of the script ($TMPDIR and also the command to create the subdirectory, below it.).

### Usage

The script accepts various command line parameters:
..* -wildseq and - -mutseq are the two oligonucleotide sequences, as strings.
..* -mutations is the number of mutations allowed to create the restriction site
..* -scan (left|right|straddle|all) controls whether just mutations to the left, right, straddling (where - -mutations > 2) or all, are reported.
..* -cut (all|wild|mut) outputs suggestions which cut either, just the wildtype or just the mutant sequence.
..* -zero includes - -mutations=0, in order to indicate confounding cut sites.
..* -revcomp reverse complements all sequences.

-h --help --perldoc give help

The script can be run at the command line, easily incorporated into a galaxy (https://usegalaxy.org/) installation, or could pretty easily be modified to run as a web cgi.

### Copyright

Copyright © 2018 a Syngenta group company

### License Information

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

THIS PROGRAM IS MADE AVAILABLE FOR DISTRIBUTION WITHOUT ANY FORM OF WARRANTY TO THE 
EXTENT PERMITTED BY APPLICABLE LAW.  THE COPYRIGHT HOLDER PROVIDES THE PROGRAM \"AS IS\" 
WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM LIES
WITH THE USER.  SHOULD THE PROGRAM PROVE DEFECTIVE IN ANY WAY, THE USER ASSUMES THE
COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION. THE COPYRIGHT HOLDER IS NOT 
RESPONSIBLE FOR ANY AMENDMENT, MODIFICATION OR OTHER ENHANCEMENT MADE TO THE PROGRAM 
BY ANY USER WHO REDISTRIBUTES THE PROGRAM SO AMENDED, MODIFIED OR ENHANCED.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL THE 
COPYRIGHT HOLDER BE LIABLE TO ANY USER FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL,
INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE
OR LOSSES SUSTAINED BY THE USER OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO 
OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER HAS BEEN ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGES.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

### Contact Details

Contact can be made at open.publishing@syngenta.com






