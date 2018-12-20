#!/usr/bin/perl -w
#
# summarize - summarize stats from *.csv files;
#             in particular, calculate out energy and area information
#
# example: summarize cornell*.csv
#
use strict;
use warnings;

@ARGV or error( "no files given" );

my $op_costs_file = "RT_FF_Analysis - CORDIC Energy and Area.csv";

my $all_stats = [];         # one entry per file
my $op_costs  = {};

for my $file ( @ARGV )
{
    read_stats( $file );
}
read_op_costs();
print_summary();


sub read_stats
{
    my $file = shift;

    my $csv = read_csv( $file );

    # skip to after "Totals:" line
    #
    my $r;
    for( $r = 0; $r < @{$csv}; $r++ )
    {
        my $col_data = $csv->[$r];
        if ( @{$col_data} != 0 && $col_data->[0] eq "Totals:" ) {
            $r++;
            last;
        }
    }

    # read each total 
    #
    my $stats = {};
    my $got_one = 0;
    for( ; $r < @{$csv}; $r++ )
    {
        my $col_data = $csv->[$r];
        @{$col_data} >= 3 or error( "bad total line" );
        my $stat = {};
        my $op = $col_data->[0];
        $stats->{$op} = $stat;
        $stat->{total} = $col_data->[1];
        $stat->{scaled_total} = $col_data->[2];
        $got_one = 1;
    }
   
    $got_one or error( "could not find totals in $file" );
    push @{$all_stats}, $stats;
}

sub read_op_costs
{
    my $csv = read_csv( $op_costs_file );

    # skip to after "Math Op" line
    #
    my $r;
    for( $r = 0; $r < @{$csv}; $r++ )
    {
        my $col_data = $csv->[$r];
        if ( @{$col_data} != 0 && $col_data->[0] eq "Math Op" ) {
            $r++;
            last;
        }
    }

}

sub read_csv
{
    my $file = shift;

    my $csv = [];
    open( C, $file ) or error( "$file could not opened for reading" );
}

sub read_csv
{
    my $file = shift;

    my $csv = [];
    open( C, $file ) or error( "$file could not opened for reading" );
    while( <C> )
    {
        # separate fields
        #
        my $line = $_;
        chomp $line;
        my $col_data = [];
        my $len = length( $line );
        my $i;
        for( $i = 0; $i < $len; )
        {
            # skip whitespace
            my $ch;
            for( ; $i < $len; $i++ )
            {
                $ch = substr( $line, $i, 1 );
                $ch ne " " && $ch ne "\t" && $ch ne "\n" and last;
            } 
            $i == $len and last;
            my $s = "";
            if ( $ch eq "\"" ) {
                # string
                $i++;
                my $got_end_delim = 0;
                while( $i != $len )
                {
                    $ch = substr( $line, $i, 1 );
                    if ( $ch eq "\"" ) {
                        $got_end_delim = 1;
                        $i++;
                        last;
                    } else {
                        if ( $ch eq "\\" ) {
                            $i++;
                            $i == $len and error( "unterminated string: $line" );
                            $ch = substr( $line, $i, 1 );
                        }
                        $s .= $ch;
                        $i++;
                    }
                }
                $got_end_delim or error( "unterminated string: $line" );

                # skip whitespace
                for( ; $i < $len; $i++ )
                {
                    $ch = substr( $line, $i, 1 );
                    $ch ne " " && $ch ne "\t" && $ch ne "\n" and last;
                } 
                $i != $len && $ch eq "," and $i++;
            } else {
                # include all chars until comma or end of line
                for( ;; )
                {
                    $s .= $ch;
                    $i++;
                    $i == $len and last;
                    $ch = substr( $line, $i, 1 );
                    defined $ch or error( "bad char at column $i: $line" );
                    if ( $ch eq "," ) {
                        $i++;
                        last;
                    }
                }
            }
            push @{$col_data}, $s;
        }

        push @{$csv}, $col_data;
    }
    close( C );
    return $csv;
}

sub print_summary
{
    @ARGV == @{$all_stats} or error( "something is wrong" );
    for my $i ( 0 .. @{$all_stats}-1 )
    {
        my $file  = $ARGV[$i];
        my $stats = $all_stats->[$i];
        print "\n$file:\n";
        for my $op ( sort keys %{$stats} ) 
        {
            my $total        = $stats->{$op}->{total};
            my $scaled_total = $stats->{$op}->{scaled_total};
            print "    $op    $total    $scaled_total\n";
        }
    }
}

sub error
{
    my $msg = shift;
    
    print "ERROR: $msg\n";
    exit( 1 );
}
