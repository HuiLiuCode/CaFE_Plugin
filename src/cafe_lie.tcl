package provide cafe_lie 1.0

namespace eval ::cafe::lie:: {
    package require cafe_tools

    namespace import ::cafe::tools::*

    # global variables for the package
    variable validargs {
        -top_bound -trj_bound -top_free -trj_free -top_type -trj_type
        -par -out -debug -first_bound -last_bound -stride_bound
        -lig_bound -first_free -last_free -stride_free -lig_free -mm_exe
        -alpha -beta -gamma
    }

    variable topfile ""
    variable trjfile { }
    variable toptype "auto"
    variable trjtype "auto"
    variable parfile { }
    variable outfile "result.log"
    variable debug 0
    variable first_bound 0
    variable last_bound -1
    variable stride_bound 1
    variable ligsel_bound ""
    variable first_free 0
    variable last_free -1
    variable stride_free 1
    variable ligsel_free ""
    variable alpha 0.18
    variable beta 0.33
    variable gamma 0.0
    variable mm_exe "namd2"
}

# print usage info
proc ::cafe::lie::print_usage { } {
    show -info "Usage: lie -top_bound filename -trj_bound filename -top_free filename -trj_free filename \[-args...\]"
    show -info "Mandatory arguments:"
    show -info "  -top_bound <topology filename for the bound state>"
    show -info "  -trj_bound <trajectory filename for the bound state>"
    show -info "  -top_free <topology filename for the free state>"
    show -info "  -trj_free <trajectory filename for the free state>"
    show -info "Optional arguments:"
    show -info "  -top_type <topology filetype>                      -- default: auto"
    show -info "  -trj_type <trajectory filetype>                    -- default: auto"
    show -info "  -par <force field parameters>                      -- default: [file join $::env(CAFEDIR) par_all22_prot.inp]"
    show -info "  -out <output filename>                             -- default: result.log"
    show -info "  -debug <debug level>                               -- default: 0"
    show -info "  -first_bound <first frame for the bound state>     -- default: 0"
    show -info "  -last_bound <last frame for the bound state>       -- default: -1"
    show -info "  -stride_bound <stride for the bound state>         -- default: 1"
    show -info "  -lig_bound <ligand selection for the bound state>  -- default: \"\""
    show -info "  -first_free <first frame for the free state>       -- default: 0"
    show -info "  -last_free <last frame for the free state>         -- default: -1"
    show -info "  -stride_free <stride for the free state>           -- default: 1"
    show -info "  -lig_free <ligand selection for the free state>    -- default: \"\""
    show -info "  -mm_exe <path to NAMD>                             -- default: \"namd2\""
    show -info "  -alpha <vdW coefficient>                           -- default: 0.18"
    show -info "  -beta <electrostatic coefficient>                  -- default: 0.33"
    show -info "  -gamma <offset>                                    -- default: 0.0"
    show ""
}

proc ::cafe::lie::lie { args } {
    variable validargs

    variable topfile_bound
    variable trjfile_bound
    variable topfile_free
    variable trjfile_free
    variable toptype
    variable trjtype
    variable parfile
    variable outfile
    variable debug
    variable first_bound
    variable last_bound
    variable stride_bound
    variable ligsel_bound
    variable first_free
    variable last_free
    variable stride_free
    variable ligsel_free

    variable mm_exe
    variable alpha
    variable beta
    variable gamma

    # *********************************************************
    # **************** Do the Preparation Work ****************
    # *********************************************************
    show -info "Sanity check"
    set start0 [clock seconds]

    # parse the command-line
    set nargs [llength $args]
    if { !$nargs || $nargs < 8 || [expr $nargs % 2] } { print_usage }

    foreach { key val } $args {
        if { [string match -?* $key] } {
            switch -nocase -- $key {
                -top_bound {
                    set topfile_bound [check_file $val "top_bound"]
                }
                -top_free {
                    set topfile_free [check_file $val "top_free"]
                }
                -top_type {
                    set toptype [check_string $val "top_type"]
                }
                -trj_bound {
                    lappend trjfile_bound [check_file $val "trj_bound"]
                }
                -trj_free {
                    lappend trjfile_free [check_file $val "trj_free"]
                }
                -trj_type {
                    set trjtype [check_string $val "trj_type"]
                }
                -par {
                    lappend parfile [check_file $val "par"]
                }
                -out {
                    set outfile $val
                }
                -debug {
                    set debug [check_int $val "debug"]
                }
                -first_bound {
                    set first_bound [check_int $val "first_bound"]
                }
                -last_bound {
                    set last_bound [check_int $val "last_bound"]
                }
                -stride_bound {
                    set stride_bound [check_pos_int $val "stride_bound"]
                }
                -lig_bound {
                    set ligsel_bound $val
                }
                -first_free {
                    set first_free [check_int $val "first_free"]
                }
                -last_free {
                    set last_free [check_int $val "last_free"]
                }
                -stride_free {
                    set stride_free [check_pos_int $val "stride_free"]
                }
                -lig_free {
                    set ligsel_free $val
                }
                -mm_exe {
                    set mm_exe $val
                }
                -alpha {
                    set alpha [check_real $val "alpha"]
                }
                -beta {
                    set beta [check_real $val "beta"]
                }
                -gamma {
                    set gamma [check_real $val "gamma"]
                }
                default {
                    if { $key ni $validargs } {
                        show -err "Found unknown argument '$key'"
                    }
                }
            }
        } else {
            show -err "Found unknown argument '$key'"
        }
    }

    # check mandatory arguments
    if { $topfile_bound eq "" || $topfile_free eq "" } {
        show -err "Need topology files for the bound and free states!"
    }

    if { ![llength $trjfile_bound] || ![llength $trjfile_free] } {
        show -err "Need trajectory files for the bound and free states!"
    }

    # check optional arguments
    if { ![llength $parfile] } {
        lappend parfile [file join $::env(CAFEDIR) "par_all22_prot.inp"]
    }

    check_exe $mm_exe "mm_exe"

    foreach name { bound free } topfile [list $topfile_bound $topfile_free] \
            selstr [list $ligsel_bound $ligsel_free] first [list $first_bound $first_free] \
            last [list $last_bound $last_free] stride [list $stride_bound $stride_free] \
            trjfile [list $trjfile_bound $trjfile_free] {
        # load topology file
        if { $toptype eq "auto" } {
            set currmol [mol new $topfile waitfor all]
        } else {
            set currmol [mol new $topfile type $toptype waitfor all]
        }

        set toptype [molinfo $currmol get filetype]

        # check topology type, currently NAMD fully supports AMBER and CHARMM/X-PLOR formats
        if { $toptype ne "psf" && $toptype ne "parm" && $toptype ne "parm7" } {
            show -err "Currently only AMBER- and CHARMM/X-PLOR-formatted topology files are supported"
        }

        # check selections
        if { $selstr eq "" } {
            show -err "Selection of ligand in the $name state should be specified"
        } else {
            set sel [atomselect $currmol $selstr]
            set natoms [$sel num]
            if { !$natoms } {
                show -err "Found zero atoms for ligand in the $name state"
            } else {
                show -info "Found $natoms atoms for ligand in the $name state"
            }
            $sel delete
        }

        # load trajectory
        if { $last != -1 && $last < $first } {
            show -err "'last' should be not less than 'first'"
        }

        foreach f $trjfile {
            if { $trjtype eq "auto" } {
                mol addfile $f waitfor all
            } else {
                mol addfile $f type $trjtype waitfor all
            }
        }

        set old_nframes [molinfo $currmol get numframes]
        show -info "Loaded $old_nframes frames for ligand in the $name state"

        # delete unnecessary frames
        if { $first > 0 } {
            animate delete end [expr $first - 1] $currmol
        }
        if { $last != -1 && $last < [expr $old_nframes - 1] } {
            animate delete beg [expr $last - $first + 1] $currmol
        }

        show -info "Generating new trajectory for ligand in the $name state"

        set tmp_trj "_lie_${name}_tmp.dcd"
        animate write dcd $tmp_trj waitfor all skip $stride $currmol

        # delete old files and reload the new trajectory to save memory
        mol delete $currmol
    }

    foreach { d h m s } [timer $start0] { break }
    show -info "It took $d days $h hrs $m min $s sec"

    # *****************************************************
    # **************** Do the Calculations ****************
    # *****************************************************
    show -info "Calculating the interaction energy"
    set start [clock seconds]

    foreach name { bound free } topfile [list $topfile_bound $topfile_free] \
            first [list $first_bound $first_free] stride [list $stride_bound $stride_free] \
            selstr [list $ligsel_bound $ligsel_free] {
        set currmol [mol new $topfile type $toptype waitfor all]

        set tmp_trj "_lie_${name}_tmp.dcd"
        mol addfile $tmp_trj type dcd waitfor all
        show -info "Loaded [molinfo $currmol get numframes] frames for ligand in the $name state"

        set mm_result [calc_mm $currmol $name $selstr $tmp_trj $first $stride $topfile]
        if { $name eq "bound" } {
            foreach { bound_ele_list bound_vdw_list } $mm_result { break }
        } else {
            foreach { free_ele_list free_vdw_list } $mm_result { break }
        }
    }

    foreach { d h m s } [timer $start] { break }
    show -info "It took $d days $h hrs $m min $s sec"

    # *****************************************************
    # **************** Generate the Result ****************
    # *****************************************************
    show -info "Generating the result"
    set start [clock seconds]

    set result [open $outfile w]

    set tfmt "   %-6s %15s %15s %15s"
    set dfmt "   %-6s %15s %15.4f"
    set sepline " [string repeat - 60]"

    puts $result $sepline
    puts $result [format $tfmt Title Frames Mean SD]
    puts $result $sepline

    write_out $result " Bound:" $sepline $bound_vdw_list $bound_ele_list
    write_out $result " Free:" $sepline $free_vdw_list $free_ele_list

    foreach { vdw_b - } [calc_stats $bound_vdw_list] { break }
    foreach { ele_b - } [calc_stats $bound_ele_list] { break }
    foreach { vdw_f - } [calc_stats $free_vdw_list] { break }
    foreach { ele_f - } [calc_stats $free_ele_list] { break }
    set vdw_d [expr $vdw_b - $vdw_f]
    set ele_d [expr $ele_b - $ele_f]

    puts $result " Delta:"
    puts $result [format $dfmt "Vdw:" "" $vdw_d]
    puts $result [format $dfmt "Elec:" "" $ele_d]
    puts $result $sepline

    set vdwsc [expr $alpha*$vdw_d]
    set elesc [expr $beta*$ele_d]
    set gbind [expr $vdwsc + $elesc + $gamma]

    puts $result " Final: (alpha=$alpha, beta=$beta, gamma=$gamma)"
    puts $result [format $dfmt "Total:" "" $gbind]
    puts $result $sepline

    puts $result " * All energy values are in kcal/mol"

    close $result

    if { $debug < 2 } {
        file delete "_lie_bound_tmp.dcd"
        file delete "_lie_free_tmp.dcd"
    }

    foreach { d h m s } [timer $start] { break }
    show -info "It took $d days $h hrs $m min $s sec"

    foreach { d h m s } [timer $start0] { break }
    show -info "Total elapsed time: $d days $h hrs $m min $s sec"
}

proc ::cafe::lie::write_out { fp title end vdw_list ele_list } {
    puts $fp $title

    write_item $fp "Vdw:" $vdw_list
    write_item $fp "Elec:" $ele_list

    puts $fp $end
}

# ######################################################################
#                            MM related
# ######################################################################
proc ::cafe::lie::calc_mm { molid prefix selstr trajname first stride topfile } {
    variable debug

    # NOTE: There is a bug in namdEnergy plugin 1.4. It was said that "skip" was
    # "number of frames to skip", so given the first frame N, the next one is
    # thought to be N + $skip + 1. Unfortunately, when the selection is not
    # updated every frame, "animate write" is used there to generate a trajectory.
    # However, "animate write" actually uses the "skip" argument the same as a
    # "stride". That is to say, if the first frame is N, the next one is actually
    # N + $skip.

    set mm_list [run_namd $molid $prefix $selstr $trajname $topfile]

    set ele_list { }
    set vdw_list { }

    foreach items $mm_list {
        foreach { - - - - ele vdw } $items { break }
        lappend ele_list $ele
        lappend vdw_list $vdw
    }

    if { $debug > 0 } {
        set tfmt "#%14s %15s %15s"
        set dfmt "%15d %15.4f %15.4f"
        set fname ${prefix}_mm.log
        set title [format $tfmt Frame Vdw Elec]
        set frames { }

        for { set i 0 } { $i < [llength $mm_list] } { incr i } {
            lappend frames [expr $first + $i*$stride]
        }

        set data { }
        lappend data $frames
        lappend data $vdw_list
        lappend data $ele_list

        write_log $fname $title $dfmt $data
    }

    return [list $ele_list $vdw_list]
}

# calculate MM energy by NAMD
proc ::cafe::lie::run_namd { molid prefix selstr trajname topfile } {
    variable mm_exe
    variable debug

    write_namd_conf $molid $prefix $selstr $trajname $topfile
    exec $mm_exe ${prefix}_mm_tmp.namd > ${prefix}_mm_tmp.log
    set result [parse_namd ${prefix}_mm_tmp.log]
    if { $debug < 2 } { cleanup ${prefix}_mm }
    return $result
}

# write a NAMD configuration file
# revised from namdenergy1.4
proc ::cafe::lie::write_namd_conf { molid prefix selstr trajname topfile } {
    variable toptype
    variable parfile

    # write the pair interaction PDB needed for NAMD
    set all [atomselect $molid all frame first]
    $all set beta 2
    set sel [atomselect $molid $selstr frame first]
    $sel set beta 1
    $sel delete
    $all writepdb ${prefix}_mm_tmp.pdb
    $all delete

    if { [catch { set namdconf [open ${prefix}_mm_tmp.namd w] } err] } {
        show -err "Failed to write NAMD config file: ${prefix}_mm_tmp.namd"
    }

    # set up topology and parameters
    if { $toptype eq "parm" || $toptype eq "parm7" } {
        puts $namdconf "amber on"
        puts $namdconf "parmfile $topfile"
        puts $namdconf "readexclusions yes"
    } else {
        puts $namdconf "structure $topfile"
        puts $namdconf "paraTypeCharmm on"
        foreach pf $parfile {
           puts $namdconf "parameters $pf"
        }
        puts $namdconf "mergeCrossterms yes"
    }

    puts $namdconf "numsteps 1"
    puts $namdconf "exclude scaled1-4"
    if { $toptype eq "parm" || $toptype eq "parm7" } {
        puts $namdconf "1-4scaling 0.8333333333333333"
        puts $namdconf "scnb 2.0"
    }
    puts $namdconf "outputname ${prefix}_mm_tmp"
    puts $namdconf "temperature 0"
    puts $namdconf "COMmotion yes"
    puts $namdconf "cutoff 999"
    puts $namdconf "dielectric 1"
    puts $namdconf "switching off"

    puts $namdconf "pairInteraction on"
    puts $namdconf "pairInteractionFile ${prefix}_mm_tmp.pdb"
    puts $namdconf "pairInteractionCol B"
    puts $namdconf "pairInteractionGroup1 1"
    puts $namdconf "pairInteractionGroup2 2"
    puts $namdconf "pairInteractionSelf off"
    puts $namdconf "coordinates ${prefix}_mm_tmp.pdb"
    puts $namdconf "set ts 0"
    puts $namdconf "coorfile open dcd $trajname"
    puts $namdconf "while \{ \!\[coorfile read\] \} \{"
    puts $namdconf "    firstTimestep \$ts"
    puts $namdconf "    run 0"
    puts $namdconf "    incr ts"
    puts $namdconf "\}"
    puts $namdconf "coorfile close"

    close $namdconf
}

