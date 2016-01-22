package provide cafe_mmpbsa 1.0

namespace eval ::cafe::mmpbsa:: {
    package require readcharmmpar ;# for ::Pararead::getvdwparam
    package require topotools ;# for "topo guessatom"
    package require cafe_tools

    namespace import ::cafe::tools::*

    # global variables for the package
    # kT/NA = RT = 1.9859 cal/K/mol * 298.15 K = 0.59210 kcal/mol
    variable kt2kc 0.592
    variable j2cal 0.239

    variable validargs {
        -top -trj -top_type -trj_type -par -out -debug -first -last -stride
        -com -rec -lig -mm -mm_exe -mm_diel -pb -pb_exe -pb_siz -pb_crg -pb_rad
        -pb_indi -pb_exdi -pb_scale -pb_perfil -pb_prbrad -pb_linit -pb_maxc
        -pb_bndcon -pb_bcfl -pb_chgm -pb_srfm -pb_swin -pb_sdens -gb
        -gb_exdi -gb_ioncon -gb_sa -gb_sagamma -sa -sa_exe -sa_rad -sa_gamma
        -sa_beta -sa_prbrad -sa_samples
    }

    variable topfile ""
    variable trjfile { }
    variable toptype "auto"
    variable trjtype "auto"
    variable parfile { }
    variable outfile "result.log"
    variable debug 0
    variable first 0
    variable last -1
    variable stride 1
    variable comsel ""
    variable recsel ""
    variable ligsel ""

    variable mm 0
    variable mm_exe "namd2"
    variable mm_diel 1.0

    variable pb 0
    variable pb_exe "delphi77"
    variable pb_siz ""
    variable pb_crg ""
    variable pb_rad "bondi"
    variable pb_indi 1.0
    variable pb_exdi 80.0
    variable pb_scale 2.0
    variable pb_perfil 80.0
    variable pb_prbrad 1.4
    variable pb_linit 1000
    variable pb_maxc 0.0001
    variable pb_bndcon 4
    variable pb_bcfl "sdh"
    variable pb_chgm "spl0"
    variable pb_srfm "smol"
    variable pb_swin 0.3
    variable pb_sdens 10.0

    variable gb 0
    variable gb_exdi 78.5
    variable gb_ioncon 0.0
    variable gb_sa 0
    variable gb_sagamma 0.005

    variable sa 0
    variable sa_exe "apbs"
    variable sa_rad "bondi"
    variable sa_gamma 0.005
    variable sa_beta 0.0
    variable sa_prbrad 1.4
    variable sa_samples 500
}

# print usage info
proc ::cafe::mmpbsa::print_usage { } {
    show -info "Usage: mmpbsa -top filename -trj filename \[-args...\]"
    show -info "Mandatory arguments:"
    show -info "  -top <topology filename>"
    show -info "  -trj <trajectory filename>"
    show -info "Optional arguments:"
    show -info "  -top_type <topology filetype>    -- default: auto"
    show -info "  -trj_type <trajectory filetype>  -- default: auto"
    show -info "  -par <force field parameters>    -- default: [file join $::env(CAFEDIR) par_all22_prot.inp]"
    show -info "  -out <output filename>           -- default: result.log"
    show -info "  -debug <debug level>             -- default: 0"
    show -info "  -first <first frame>             -- default: 0"
    show -info "  -last <last frame>               -- default: -1"
    show -info "  -stride <stride>                 -- default: 1"
    show -info "  -com <complex selection>         -- default: \"\""
    show -info "  -rec <receptor selection>        -- default: \"\""
    show -info "  -lig <ligand selection>          -- default: \"\""
    show -info "  -mm <do gas-phase calculation>   -- default: 0"
    show -info "  -mm_exe <path to NAMD>           -- default: \"namd2\""
    show -info "  -mm_diel <dielectric constant>   -- default: 1.0"
    show -info "  -pb <do PB calculation>          -- default: 0"
    show -info "  -pb_exe <path to DelPhi/APBS>    -- default: \"delphi77\""
    show -info "  -pb_siz <radii parameter file>   -- default: \"\""
    show -info "  -pb_crg <charge parameter file>  -- default: \"\""
    show -info "  -pb_rad <type of PB radii>       -- default: bondi"
    show -info "  -pb_indi <internal dielectric>   -- default: 1.0"
    show -info "  -pb_exdi <external dielectric>   -- default: 80.0"
    show -info "  -pb_scale <scale>                -- default: 2.0"
    show -info "  -pb_perfil <percentage of fill>  -- default: 80.0"
    show -info "  -pb_prbrad <radius of probe>     -- default: 1.4"
    show -info "  -pb_linit <linear iterations>    -- default: 1000"
    show -info "  -pb_maxc <convergence threshold> -- default: 0.0001"
    show -info "  -pb_bndcon <boundary condition>  -- default: 4"
    show -info "  -pb_bcfl <boundary condition>    -- default: sdh"
    show -info "  -pb_chgm <charge method>         -- default: spl0"
    show -info "  -pb_srfm <surface method>        -- default: smol"
    show -info "  -pb_swin <spline window width>   -- default: 0.3"
    show -info "  -pb_sdens <number of grids>      -- default: 10.0"
    show -info "  -gb <do GB calculation>          -- default: 0"
    show -info "  -gb_exdi <external dielectric>   -- default: 78.5"
    show -info "  -gb_ioncon <ion concentration>   -- default: 0.0"
    show -info "  -gb_sa <do LCPO calculation>     -- default: 0"
    show -info "  -gb_sagamma <surface tension>    -- default: 0.005"
    show -info "  -sa <do SA calculation>          -- default: 0"
    show -info "  -sa_exe <path to APBS>           -- default: \"apbs\""
    show -info "  -sa_rad <type of SA radii>       -- default: bondi"
    show -info "  -sa_gamma <surface tension>      -- default: 0.005"
    show -info "  -sa_beta <surface offset>        -- default: 0.0"
    show -info "  -sa_prbrad <radius of probe>     -- default: 1.4"
    show -info "  -sa_samples <number of samples>  -- default: 500"
    show ""
}

proc ::cafe::mmpbsa::mmpbsa { args } {
    variable validargs

    variable topfile
    variable trjfile
    variable toptype
    variable trjtype
    variable parfile
    variable outfile
    variable debug
    variable first
    variable last
    variable stride
    variable comsel
    variable recsel
    variable ligsel

    variable mm
    variable mm_exe
    variable mm_diel

    variable pb
    variable pb_exe
    variable pb_siz
    variable pb_crg
    variable pb_rad
    variable pb_indi
    variable pb_exdi
    variable pb_scale
    variable pb_perfil
    variable pb_prbrad
    variable pb_linit
    variable pb_maxc
    variable pb_bndcon
    variable pb_bcfl
    variable pb_chgm
    variable pb_srfm
    variable pb_swin
    variable pb_sdens

    variable gb
    variable gb_exdi
    variable gb_ioncon
    variable gb_sa
    variable gb_sagamma

    variable sa
    variable sa_exe
    variable sa_rad
    variable sa_gamma
    variable sa_beta
    variable sa_prbrad
    variable sa_samples

    # *********************************************************
    # **************** Do the Preparation Work ****************
    # *********************************************************
    show -info "Sanity check"
    set start0 [clock seconds]

    # parse the command-line
    set nargs [llength $args]
    if { !$nargs || $nargs < 4 || [expr $nargs % 2] } { print_usage }

    foreach { key val } $args {
        if { [string match -?* $key] } {
            switch -nocase -- $key {
                -top {
                    set topfile [check_file $val "top"]
                }
                -top_type {
                    set toptype [check_string $val "top_type"]
                }
                -trj {
                    lappend trjfile [check_file $val "trj"]
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
                -first {
                    set first [check_int $val "first"]
                }
                -last {
                    set last [check_int $val "last"]
                }
                -stride {
                    set stride [check_pos_int $val "stride"]
                }
                -com {
                    set comsel $val
                }
                -rec {
                    set recsel $val
                }
                -lig {
                    set ligsel $val
                }
                -mm {
                    set mm [check_bool $val "mm"]
                }
                -pb {
                    set pb [check_nneg_int $val "pb"]
                }
                -gb {
                    set gb [check_nneg_int $val "gb"]
                }
                -sa {
                    set sa [check_nneg_int $val "sa"]
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

    if { !$mm && !$pb && !$gb && !$sa } {
        show -err "Need a calculation type!"
    }

    if { $mm } {
        foreach { key val } $args {
            switch -nocase -- $key {
                -mm_exe {
                    set mm_exe $val
                }
                -mm_diel {
                    set mm_diel [check_pos_real $val "mm_diel"]
                }
                default {
                    continue
                }
            }
        }
    }

    if { $pb } {
        foreach { key val } $args {
            switch -nocase -- $key {
                -pb_exe {
                    set pb_exe $val
                }
                -pb_siz {
                    set pb_siz [check_file $val "pb_siz"]
                }
                -pb_crg {
                    set pb_crg [check_file $val "pb_crg"]
                }
                -pb_rad {
                    set pb_rad [check_string $val "pb_rad"]
                }
                -pb_indi {
                    set pb_indi [check_pos_real $val "pb_indi"]
                }
                -pb_exdi {
                    set pb_exdi [check_pos_real $val "pb_exdi"]
                }
                -pb_scale {
                    set pb_scale [check_pos_real $val "pb_scale"]
                }
                -pb_perfil {
                    set pb_perfil [check_pos_real $val "pb_perfil"]
                    if { $pb_perfil > 100 } {
                        show -err "'pb_perfil' should be less than 100"
                    }
                }
                -pb_prbrad {
                    set pb_prbrad [check_nneg_real $val "pb_prbrad"]
                }
                -pb_linit {
                    set pb_linit [check_pos_int $val "pb_linit"]
                }
                -pb_maxc {
                    set pb_maxc [check_pos_real $val "pb_maxc"]
                }
                -pb_bndcon {
                    set pb_bndcon [check_pos_int $val "pb_bndcon"]
                }
                -pb_bcfl {
                    set pb_bcfl [check_string $val "pb_bcfl"]
                }
                -pb_chgm {
                    set pb_chgm [check_string $val "pb_chgm"]
                }
                -pb_srfm {
                    set pb_srfm [check_string $val "pb_srfm"]
                }
                -pb_swin {
                    set pb_swin [check_pos_real $val "pb_swin"]
                }
                -pb_sdens {
                    set pb_sdens [check_pos_real $val "pb_sdens"]
                }
                default {
                    continue
                }
            }
        }
    }

    if { $gb } {
        foreach { key val } $args {
            switch -nocase -- $key {
                -mm_exe {
                    set mm_exe $val
                }
                -gb_exdi {
                    set gb_exdi [check_pos_real $val "gb_exdi"]
                }
                -gb_ioncon {
                    set gb_ioncon [check_nneg_real $val "gb_ioncon"]
                }
                -gb_sa {
                    set gb_sa [check_bool $val "gb_sa"]
                }
                -gb_sagamma {
                    set gb_sagamma [check_real $val "gb_sagamma"]
                }
                default {
                    continue
                }
            }
        }
    }

    if { $sa } {
        foreach { key val } $args {
            switch -nocase -- $key {
                -sa_exe {
                    set sa_exe $val
                }
                -sa_rad {
                    set sa_rad [check_string $val "sa_rad"]
                }
                -sa_gamma {
                    set sa_gamma [check_real $val "sa_gamma"]
                }
                -sa_beta {
                    set sa_beta [check_real $val "sa_beta"]
                }
                -sa_prbrad {
                    set sa_prbrad [check_pos_real $val "sa_prbrad"]
                }
                -sa_samples {
                    set sa_samples [check_pos_int $val "sa_samples"]
                }
                default {
                    continue
                }
            }
        }
    }

    # check mandatory arguments
    if { $topfile eq "" } {
        show -err "Need a topology file!"
    }

    if { ![llength $trjfile] } {
        show -err "Need a trajectory file!"
    }

    # check optional arguments
    if { ![llength $parfile] } {
        lappend parfile [file join $::env(CAFEDIR) "par_all22_prot.inp"]
    }

    if { $mm } { check_exe $mm_exe "mm_exe"}
    if { $pb } { check_exe $pb_exe "pb_exe"}
    if { $gb } { check_exe $mm_exe "mm_exe"}
    if { $sa == 2 } { check_exe $sa_exe "sa_exe"}

    if { $pb && $gb } {
        show -err "Select only one implicit model at a time."
    }

    if { $gb && ($gb_sa && $sa) } {
        show -err "Select only one SASA method at a time."
    }

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
    if { $comsel eq "" } {
        show -err "Selection of complex must be specified"
    }

    if { ($recsel eq "" && $ligsel ne "") || ($recsel ne "" && $ligsel eq "") } {
        show -err "Selections of both receptor and ligand should be specified"
    }

    foreach name { complex receptor ligand } selstr [list $comsel $recsel $ligsel] {
        if { $selstr ne "" } {
            set sel [atomselect $currmol $selstr]
            set natoms [$sel num]
            if { !$natoms } {
                show -err "Found zero atoms for $name"
            } else {
                show -info "Found $natoms atoms for $name"
            }
            $sel delete
        }
    }

    # check radii setting
    if { $pb } {
        if { $pb_rad eq "charmm" && ![llength $parfile] } {
            show -err "Need a CHARMM formatted parameter file for pb_rad charmm"
        } elseif { $pb_rad eq "parm7" && ($toptype ne "parm" && $toptype ne "parm7") } {
            show -err "Need an AMBER PARM7 file for pb_rad parm7"
        }
    }

    if { $sa } {
        if { $sa_rad eq "charmm" && ![llength $parfile] } {
            show -err "Need a CHARMM formatted parameter file for sa_rad charmm"
        } elseif { $sa_rad eq "parm7" && ($toptype ne "parm" && $toptype ne "parm7") } {
            show -err "Need an AMBER PARM7 file for sa_rad parm7"
        }
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
    show -info "Loaded $old_nframes frames for complex"

    # delete unnecessary frames
    if { $first > 0 } {
        animate delete end [expr $first - 1] $currmol
    }
    if { $last != -1 && $last < [expr $old_nframes - 1] } {
        animate delete beg [expr $last - $first + 1] $currmol
    }

    show -info "Generating new trajectory for complex"

    set com_trj "_mmpbsa_com_tmp.dcd"

    animate write dcd $com_trj waitfor all skip $stride $currmol

    # delete old files and reload the new trajectory to save memory
    mol delete $currmol

    set currmol [mol new $topfile type $toptype waitfor all]

    mol addfile $com_trj type dcd waitfor all
    show -info "Loaded [molinfo $currmol get numframes] frames for complex"

    foreach { d h m s } [timer $start0] { break }
    show -info "It took $d days $h hrs $m min $s sec"

    # ***************************************************************
    # **************** Do the Gas-Phase Calculations ****************
    # ***************************************************************
    if { $mm } {
        show -info "Calculating the MM term"
        set start [clock seconds]

        set com_mm_result [calc_mm $currmol com $comsel $com_trj]
        foreach { com_ele_list com_vdw_list } $com_mm_result { break }

        if { $recsel ne "" } {
            set rec_mm_result [calc_mm $currmol rec $recsel $com_trj]
            foreach { rec_ele_list rec_vdw_list } $rec_mm_result { break }
        }

        if { $ligsel ne "" } {
            set lig_mm_result [calc_mm $currmol lig $ligsel $com_trj]
            foreach { lig_ele_list lig_vdw_list } $lig_mm_result { break }
        }

        foreach { d h m s } [timer $start] { break }
        show -info "It took $d days $h hrs $m min $s sec"
    }

    # *********************************************************************
    # **************** Do the Polar Solvation Calculations ****************
    # *********************************************************************
    if { $pb } {
        show -info "Calculating the PB term"
        set start [clock seconds]

        # assign radii
        set ar_args "$currmol $pb_rad"
        if { $pb_rad eq "charmm" } {
            foreach p $parfile { append ar_args " \"$p\"" }
        } elseif { $pb_rad eq "parm7" } {
            append ar_args " \"$topfile\""
        }
        eval assign_radii $ar_args

        set com_pb_list [calc_pb $currmol com $comsel]

        if { $recsel ne "" } {
            set rec_pb_list [calc_pb $currmol rec $recsel]
        }

        if { $ligsel ne "" } {
            set lig_pb_list [calc_pb $currmol lig $ligsel]
        }

        foreach { d h m s } [timer $start] { break }
        show -info "It took $d days $h hrs $m min $s sec"
    }

    # ************************************************************************
    # **************** Do the Nonpolar Solvation Calculations ****************
    # ************************************************************************
    if { $sa } {
        show -info "Calculating the SA term"
        set start [clock seconds]

        # assign radii again
        set ar_args "$currmol $sa_rad"
        if { $sa_rad eq "charmm" } {
            foreach p $parfile { append ar_args " \"$p\"" }
        } elseif { $sa_rad eq "parm7" } {
            append ar_args " \"$topfile\""
        }
        eval assign_radii $ar_args

        set com_sa_list [calc_sa $currmol com $comsel]

        if { $recsel ne "" } {
            set rec_sa_list [calc_sa $currmol rec $recsel]
        }

        if { $ligsel ne "" } {
            set lig_sa_list [calc_sa $currmol lig $ligsel]
        }

        foreach { d h m s } [timer $start] { break }
        show -info "It took $d days $h hrs $m min $s sec"
    }

    # *****************************************************
    # **************** Generate the Result ****************
    # *****************************************************
    show -info "Generating result"
    set start [clock seconds]

    set result [open $outfile w]

    set tfmt "   %-6s %15s %15s %15s"
    set sepline " [string repeat - 60]"

    puts $result $sepline
    puts $result [format $tfmt Title Frames Mean SD]
    puts $result $sepline

    if { !$mm } {
        set com_ele_list { }
        set com_vdw_list { }
    }
    if { !$pb } {
        set com_pb_list { }
    }
    if { !$sa } {
        set com_sa_list { }
    }

    write_out $result " Complex:" $sepline $com_ele_list $com_vdw_list $com_pb_list $com_sa_list

    if { $recsel ne "" } {
        if { !$mm } {
            set rec_ele_list { }
            set rec_vdw_list { }
        }
        if { !$pb } {
            set rec_pb_list { }
        }
        if { !$sa } {
            set rec_sa_list { }
        }

        write_out $result " Receptor:" $sepline $rec_ele_list $rec_vdw_list $rec_pb_list $rec_sa_list
    }

    if { $ligsel ne "" } {
        if { !$mm } {
            set lig_ele_list { }
            set lig_vdw_list { }
        }
        if { !$pb } {
            set lig_pb_list { }
        }
        if { !$sa } {
            set lig_sa_list { }
        }

        write_out $result " Ligand:" $sepline $lig_ele_list $lig_vdw_list $lig_pb_list $lig_sa_list
    }

    if { $comsel ne "" && $recsel ne "" && $ligsel ne "" } {
        set del_ele_list [vecsub $com_ele_list [vecadd $rec_ele_list $lig_ele_list]]

        set del_vdw_list [vecsub $com_vdw_list [vecadd $rec_vdw_list $lig_vdw_list]]

        set del_pb_list [vecsub $com_pb_list [vecadd $rec_pb_list $lig_pb_list]]

        set del_sa_list [vecsub $com_sa_list [vecadd $rec_sa_list $lig_sa_list]]

        write_out $result " Delta:" $sepline $del_ele_list $del_vdw_list $del_pb_list $del_sa_list
    }

    puts $result " * All energy values are in kcal/mol"

    close $result

    if { $debug < 2 } { file delete $com_trj }

    foreach { d h m s } [timer $start] { break }
    show -info "It took $d days $h hrs $m min $s sec"

    foreach { d h m s } [timer $start0] { break }
    show -info "Total elapsed time: $d days $h hrs $m min $s sec"
}

proc ::cafe::mmpbsa::write_out { fp title end ele_list vdw_list pb_list sa_list } {
    variable mm
    variable pb
    variable sa

    puts $fp $title

    if { $mm } {
        write_item $fp "Elec:" $ele_list
        write_item $fp "Vdw:" $vdw_list
    }

    if { $pb } {
        write_item $fp "PB:" $pb_list
    }

    if { $sa } {
        write_item $fp "SA:" $sa_list
    }

    if { $mm } {
        set gas_list [vecadd $ele_list $vdw_list]
        write_item $fp "Gas:" $gas_list
    }

    if { $pb && $sa } {
        set sol_list [vecadd $pb_list $sa_list]
        write_item $fp "Sol:" $sol_list
    }

    if { $mm && $pb } {
        set pol_list [vecadd $ele_list $pb_list]
        write_item $fp "Pol:" $pol_list
    }

    if { $mm && $sa } {
        set npol_list [vecadd $vdw_list $sa_list]
        write_item $fp "Npol:" $npol_list
    }

    set tot_list { }
    if { $mm } {
        set tot_list $gas_list
    }
    if { $pb } {
        if { ![llength $tot_list] } {
            set tot_list $pb_list
        } else {
            set tot_list [vecadd $tot_list $pb_list]
        }
    }
    if { $sa } {
        if { ![llength $tot_list] } {
            set tot_list $sa_list
        } else {
            set tot_list [vecadd $tot_list $sa_list]
        }
    }
    write_item $fp "Total:" $tot_list
    puts $fp $end
}

# ######################################################################
#                            MM related
# ######################################################################
proc ::cafe::mmpbsa::calc_mm { molid prefix selstr trajname } {
    variable first
    variable stride
    variable debug

    # NOTE: There is a bug in namdEnergy plugin 1.4. It was said that "skip" was
    # "number of frames to skip", so given the first frame N, the next one is
    # thought to be N + $skip + 1. Unfortunately, when the selection is not
    # updated every frame, "animate write" is used there to generate a trajectory.
    # However, "animate write" actually uses the "skip" argument the same as a
    # "stride". That is to say, if the first frame is N, the next one is actually
    # N + $skip.

    set mm_list [run_namd $molid $prefix $selstr $trajname]

    set ele_list { }
    set vdw_list { }
    set tot_list { }

    foreach items $mm_list {
        foreach { - - - - ele vdw } $items { break }
        lappend ele_list $ele
        lappend vdw_list $vdw
        lappend tot_list [expr $ele + $vdw]
    }

    if { $debug > 0 } {
        set tfmt "#%14s %15s %15s %15s"
        set dfmt "%15d %15.4f %15.4f %15.4f"
        set fname ${prefix}_mm.log
        set title [format $tfmt Frame Elec Vdw MM]
        set frames { }

        for { set i 0 } { $i < [llength $mm_list] } { incr i } {
            lappend frames [expr $first + $i*$stride]
        }

        set data { }
        lappend data $frames
        lappend data $ele_list
        lappend data $vdw_list
        lappend data $tot_list

        write_log $fname $title $dfmt $data
    }

    return [list $ele_list $vdw_list]
}

# calculate MM and/or GB by NAMD
proc ::cafe::mmpbsa::run_namd { molid prefix selstr trajname } {
    variable mm_exe
    variable debug

    write_namd_conf $molid $prefix $selstr $trajname
    exec $mm_exe ${prefix}_mm_tmp.namd > ${prefix}_mm_tmp.log
    set result [parse_namd ${prefix}_mm_tmp.log]
    if { $debug < 2 } { cleanup ${prefix}_mm }
    return $result
}

# write a namd configuration file
# revised from namdenergy1.4
proc ::cafe::mmpbsa::write_namd_conf { molid prefix selstr trajname } {
    variable toptype
    variable topfile
    variable parfile
    variable mm_diel
    variable gb
    variable gb_exdi
    variable gb_ioncon
    variable gb_sa
    variable gb_sagamma

    # write the pair interaction PDB needed for NAMD
    set all [atomselect $molid all frame first]
    $all set beta 0
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
    puts $namdconf "dielectric $mm_diel"
    puts $namdconf "switching off"

    if { $gb } {
        puts $namdconf "GBIS on"
        puts $namdconf "solventDielectric $gb_exdi"
        puts $namdconf "intrinsicRadiusOffset 0.09"
        puts $namdconf "ionConcentration $gb_ioncon"
        puts $namdconf "alphaCutoff 999"
        if { $gb == 1 } {
            puts $namdconf "GBISDelta 1.0"
            puts $namdconf "GBISBeta 0.8"
            puts $namdconf "GBISGamma 4.85"
        } elseif { $gb == 2 } {
            puts $namdconf "GBISDelta 0.8"
            puts $namdconf "GBISBeta 0.0"
            puts $namdconf "GBISGamma 2.90912"
        } else {
            show -err "Unknown GB method: $gb"
        }
        if { $gb_sa } {
            puts $namdconf "SASA on"
            puts $namdconf "surfaceTension $gb_sagamma"
        } else {
            puts $namdconf "SASA off"
        }
    } else {
        puts $namdconf "GBIS off"
    }

    puts $namdconf "pairInteraction on"
    puts $namdconf "pairInteractionFile ${prefix}_mm_tmp.pdb"
    puts $namdconf "pairInteractionCol B"
    puts $namdconf "pairInteractionGroup1 1"
    puts $namdconf "pairInteractionSelf on"
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

# ######################################################################
#                              PB related
# ######################################################################
proc ::cafe::mmpbsa::assign_radii { molid type args } {
    switch $type {
        bondi     { assign_radii_bondi $molid }
        rowland   { assign_radii_rowland $molid }
        mparse    { assign_radii_parse $molid -mod }
        parse     { assign_radii_parse $molid -orig }
        charmm    { eval assign_radii_charmm $molid $args }
        roux      { assign_radii_roux $molid }
        parm7     { assign_radii_parm7 $molid $args }
        swanson   { assign_radii_swanson $molid }
        tan       { assign_radii_tan $molid }
        yamagishi { assign_radii_yamagishi $molid }
        default   { show -err "Unknown type of radii: $type" }
    }
}

# **************** General ****************

# assign Bondi radii (J Phys Chem 1964, 68, 441-452)
proc ::cafe::mmpbsa::assign_radii_bondi { molid } {
    # element guessed from mass should be more accurate; however,
    # there is a bug in topotools1.5 plugin, line 43:
    # "$s set element [lindex $elements [ptefrommass $idx]]"
    # should be
    # "$s set element [lindex $elements [ptefrommass $a]]"
    # so I have to guess from name for now ...
    topo guessatom element mass -molid $molid
    #topo guessatom element name -molid $molid
    topo guessatom radius element -molid $molid
}

# Rowland and Taylor (J Phys Chem 1996, 100, 7384-7391)
proc ::cafe::mmpbsa::assign_radii_rowland { molid } {
    assign_radii_bondi $molid

    [atomselect $molid "element 'H'"] set radius 1.10
    [atomselect $molid "element 'C'"] set radius 1.77
    [atomselect $molid "element 'N'"] set radius 1.64
    [atomselect $molid "element 'O'"] set radius 1.58
    [atomselect $molid "element 'F'"] set radius 1.46
    [atomselect $molid "element 'S'"] set radius 1.81
    [atomselect $molid "element 'Cl'"] set radius 1.76
    [atomselect $molid "element 'Br'"] set radius 1.87
    [atomselect $molid "element 'I'"] set radius 2.03
}

# assign PARSE radii (J Phys Chem 1994, 98, 1978-1988)
# atoms not specified here will use default values in VMD
proc ::cafe::mmpbsa::assign_radii_parse { molid { opt "-mod" } } {
    # unknown atoms were set to 1.5 in AMBER; however, I set them
    # according to either Pauling or the value guessed by VMD
    #[atomselect $molid "all"] set radius 1.50
    assign_radii_bondi $molid

    [atomselect $molid "element 'H'"] set radius 1.00
    [atomselect $molid "element 'C'"] set radius 1.70
    [atomselect $molid "element 'O'"] set radius 1.40
    [atomselect $molid "element 'N'"] set radius 1.50
    [atomselect $molid "element 'S'"] set radius 1.85

    # J Mol Biol 2007, 366, 1475-1496
    [atomselect $molid "element 'P'"] set radius 1.9

    # Pauling radii, from which PARSE is derived
    [atomselect $molid "element 'F'"] set radius 1.35
    [atomselect $molid "element 'Cl'"] set radius 1.80
    [atomselect $molid "element 'Br'"] set radius 1.95
    [atomselect $molid "element 'I'"] set radius 2.15

    # set CHn = 2.0 and nonpolar H = 0.0
    if { $opt == "-orig" } {
        return
        # the application is wrong for aromatic ring
        set hs [atomselect $molid "element 'H' and numbonds = 1"]
        set nphlist { }
        set chnlist { }
        set i 0

        foreach bid [$hs getbonds] {
            set batom [atomselect $molid "index $bid"]
            if { [$batom get element] eq "C" } {
                lappend nphlist $i
                lappend chnlist $bid
            }
            $batom delete
            incr i
        }

        set ids [$hs get index]
        $hs delete
        foreach i $nphlist {
            [atomselect $molid "index [lindex $ids $i]"] set radius 0.00
        }

        foreach i [lsort -unique -integer $chnlist] {
            [atomselect $molid "index $i"] set radius 2.00
        }
    }
}

# **************** For CHARMM ****************

# assign CHARMM radii based on CHARMM force field parameters
proc ::cafe::mmpbsa::assign_radii_charmm { molid args } {
    set paramlist { }

    foreach f $args {
        lappend paramlist [::Pararead::read_charmm_parameters $f]
    }

    set sel [atomselect $molid all]
    set typelist [lsort -uniq [$sel get type]]
    $sel delete

    foreach type $typelist {
        set typesel [atomselect $molid "type $type"]

        if { [$typesel num] > 0 } {
            foreach par $paramlist {
                set vdwparams [::Pararead::getvdwparam $par $type]

                if { [llength $vdwparams] > 1 } {
                    set r [expr [lindex $vdwparams 1] * 1.0]
                    $typesel set radius $r
                    break;
                }
            }
        }
        $typesel delete
    }
}

# assign modified CHARMM radii developed by Roux and coworkers
# (J Phys Chem B 1997, 101, 5239-5248 & J Phys Chem B 2002, 106, 11026-11035)
proc ::cafe::mmpbsa::assign_radii_roux { molid } {
    # **************** General ****************
    assign_radii_bondi $molid

    # zero H and set heavy atoms to average default values
    [atomselect $molid "hydrogen"] set radius 0.0
    [atomselect $molid "carbon"] set radius 2.3
    [atomselect $molid "oxygen"] set radius 1.8
    [atomselect $molid "nitrogen"] set radius 2.3
    [atomselect $molid "sulfur"] set radius 2.3

    # **************** Protein ****************

    # Patches
    # CT3 N-Methylamide C-terminus
    # ACE acetylated N-terminus (ACP for PRO)
    [atomselect $molid "name CAY CAT"] set radius 2.06
    [atomselect $molid "name CY"] set radius 2.04
    [atomselect $molid "name OY"] set radius 1.52
    [atomselect $molid "name NT"] set radius 2.23
    [atomselect $molid "name \"OT.*\""] set radius 1.40 ;# for COO- terminus

    # Backbone
    [atomselect $molid "name C"] set radius 2.04 ;# for peptide bond
    [atomselect $molid "name O"] set radius 1.52 ;# for peptide bond
    [atomselect $molid "name N"] set radius 2.23 ;# for peptide bond
    [atomselect $molid "name CA"] set radius 2.86 ;# for all CA except GLY
    [atomselect $molid "resname GLY and name CA"] set radius 2.38 ;# for GLY only

    # C
    [atomselect $molid "name CB"] set radius 2.67 ;# for all residues
    [atomselect $molid "name \"CG.*\""] set radius 2.46 ;# for ARG, GLN, ILE, LYS, MET, PHE, THR, TRP, VAL, HSP, HSD
    [atomselect $molid "resname GLU and name CG"] set radius 2.77 ;# for GLU only
    [atomselect $molid "name \"CD.*\""] set radius 2.44 ;# for ARG, ILE, LEU, LYS
    [atomselect $molid "(resname GLN and name CD) or (resname ASN and name CG) or (resname GLU and name CD) or (resname ASP and name CG)"] set radius 1.98 ;# for ASP, GLU, ASN, GLN
    [atomselect $molid "resname PRO and (name CB CG CD)"] set radius 1.98 ;# for PRO only
    [atomselect $molid "(resname TYR and (name \"CE.*\" \"CD.*\" CZ)) or (resname PHE and (name \"CE.*\" \"CD.*\" CZ))"] set radius 2.00 ;# for TYR, PHE rings
    [atomselect $molid "resname TRP and (name \"CE.*\" \"CD.*\" \"CZ.*\" CH2)"] set radius 1.78 ;# for TRP ring only
    [atomselect $molid "name CE"] set radius 2.10 ;# for MET only
    [atomselect $molid "(resname ARG and name CZ) or (resname LYS and name CE)"] set radius 2.80 ;# for ARG, LYS
    #[atomselect $molid "((resname HSD HSP) and name CE1) or ((resname HSD HSP) and name CD2)"] set radius 1.98 ;# for neutral HSD and protonated HSP
    [atomselect $molid "((resname HSD HSP HSE HIS) and name CE1) or ((resname HSD HSP HSE HIS) and name CD2)"] set radius 1.98 ;# for neutral HSD and protonated HSP, and HSE

    # O
    [atomselect $molid "(resname GLU ASP) and (name \"OE.*\" \"OD.*\")"] set radius 1.40 ;# for GLU, ASP
    [atomselect $molid "(resname ASN GLN) and (name \"OE.*\" \"OD.*\")"] set radius 1.42 ;# for ASN, GLN
    [atomselect $molid "name \"OG.*\""] set radius 1.64 ;# for SER, THR
    [atomselect $molid "resname TYR and name OH"] set radius 1.85 ;# for TYR only
    [atomselect $molid "resname TIP3 and name OH2"] set radius 2.2 ;# for explicit water molecules

    # N
    #[atomselect $molid "resname HSD and (name NE2 ND1)"] set radius 1.80 ;# for neutral HSD
    [atomselect $molid "(resname HSD HSE HIS) and (name NE2 ND1)"] set radius 1.80 ;# for neutral HSD, and HSE
    [atomselect $molid "resname HSP and (name NE2 ND1)"] set radius 2.30 ;# for protonated HSP
    [atomselect $molid "(resname ARG and (name \"NH.*\" NE)) or (resname LYS and name NZ)"] set radius 2.13 ;# for ARG, LYS
    [atomselect $molid "(resname GLN and name NE2) or (resname ASN and name ND2)"] set radius 2.15 ;# for GLN, ASN
    [atomselect $molid "resname TRP and name NE1"] set radius 2.40 ;# for TRP

    # S
    [atomselect $molid "name \"S.*\" and (resname MET CYS)"] set radius 2.00 ;# for MET, CYS

    # **************** Ions ****************

    [atomselect $molid "resname POT"] set radius 2.035 ;# potassium ion K+
    [atomselect $molid "resname CLA"] set radius 2.035 ;# chloride ion Cl-
    [atomselect $molid "resname SOD"] set radius 1.66 ;# sodium ion Na+

    # Tetramethylamonium (TEA)
    [atomselect $molid "resname TEA and name N"] set radius 2.15
    [atomselect $molid "resname TEA and (name C1 C2 C3 C4)"] set radius 2.30
    [atomselect $molid "resname TEA and (name C5 C6 C7 C8)"] set radius 2.30

    # **************** Nucleic acids ****************

    # purine base atoms
    [atomselect $molid "resname ADE and name N1"] set radius 1.75
    [atomselect $molid "resname GUA and name N1"] set radius 2.15
    [atomselect $molid "(resname GUA ADE) and name C2"] set radius 2.15
    [atomselect $molid "(resname GUA ADE) and name N3"] set radius 1.69
    [atomselect $molid "(resname GUA ADE) and name C4"] set radius 2.12
    [atomselect $molid "(resname GUA ADE) and name C5"] set radius 2.12
    [atomselect $molid "(resname GUA ADE) and name C6"] set radius 2.12
    [atomselect $molid "(resname GUA ADE) and name N7"] set radius 1.69
    [atomselect $molid "(resname GUA ADE) and name C8"] set radius 2.12
    [atomselect $molid "(resname GUA ADE) and name N9"] set radius 2.13
    [atomselect $molid "resname ADE and name N6"] set radius 2.17
    [atomselect $molid "resname GUA and name O6"] set radius 1.55
    [atomselect $molid "resname GUA and name N2"] set radius 2.12
    [atomselect $molid "(resname GUA ADE) and name C9"] set radius 2.30

    # pyrimidine base atoms
    [atomselect $molid "(resname CYT THY URA) and name N1"] set radius 2.20
    [atomselect $molid "(resname CYT THY URA) and name C2"] set radius 2.04
    [atomselect $molid "resname CYT and name N3"] set radius 1.68
    [atomselect $molid "(resname THY URA) and name N3"] set radius 2.20
    [atomselect $molid "(resname CYT THY URA) and name C4"] set radius 2.12
    [atomselect $molid "(resname CYT THY URA) and name C5"] set radius 2.25
    [atomselect $molid "(resname CYT THY URA) and name C6"] set radius 2.25
    [atomselect $molid "(resname CYT THY URA) and name O2"] set radius 1.60
    [atomselect $molid "resname CYT and name N4"] set radius 2.08
    [atomselect $molid "(resname THY URA) and name O4"] set radius 1.60
    [atomselect $molid "resname THY and name C5M"] set radius 2.30
    [atomselect $molid "(resname CYT THY URA) and name C1"] set radius 2.30

    # sugar atoms
    [atomselect $molid "name C1'"] set radius 2.57
    [atomselect $molid "name C2'"] set radius 2.70
    [atomselect $molid "name C3'"] set radius 2.73
    [atomselect $molid "name C4'"] set radius 2.50
    [atomselect $molid "name O4'"] set radius 1.55
    [atomselect $molid "name C5'"] set radius 2.57
    [atomselect $molid "name O3'"] set radius 1.65
    [atomselect $molid "name O5'"] set radius 1.65
    [atomselect $molid "name O2'"] set radius 1.75

    # add radii for blocking group hydroxyl oxygens
    set oter [atomselect $molid "(name O3' O5') and numbonds = 2"]
    set oterlist { }
    set i 0

    foreach bids [$oter getbonds] {
        set batom [atomselect $molid "index [lindex bids 0] or index [lindex bids 1]"]
        set bname [$batom get name]
        if { "H3T" in $bname || "H5T" in $bname } { lappend oterlist $i }
        $batom delete
        incr i
    }

    set ids [$oter get index]
    $oter delete

    foreach i $oterlist {
        [atomselect $molid "index [lindex $ids $i]"] set radius 1.72
    }

    # atoms for sugar2phos
    [atomselect $molid "resname ADE and name O5T"] set radius 1.65
    [atomselect $molid "resname ADE and name C5T"] set radius 2.30
    [atomselect $molid "resname ADE and name P3"] set radius 2.35
    [atomselect $molid "resname ADE and name O1P3"] set radius 1.49
    [atomselect $molid "resname ADE and name O2P3"] set radius 1.49
    [atomselect $molid "resname ADE and name O3T"] set radius 1.65
    [atomselect $molid "resname ADE and name C3T"] set radius 2.30

    # phosphate atoms
    [atomselect $molid "(resname GUA ADE CYT THY URA) and name P"] set radius 2.35
    [atomselect $molid "(resname GUA ADE CYT THY URA) and name O1P"] set radius 1.49
    [atomselect $molid "(resname GUA ADE CYT THY URA) and name O2P"] set radius 1.49

    # DMPA phosphate atoms, from above radii
    [atomselect $molid "resname DMPA and name P1"] set radius 2.35
    [atomselect $molid "resname DMPA and name O1"] set radius 1.65
    [atomselect $molid "resname DMPA and name O2"] set radius 1.65
    [atomselect $molid "resname DMPA and name O3"] set radius 1.49
    [atomselect $molid "resname DMPA and name O4"] set radius 1.49
    [atomselect $molid "resname DMPA and name C1"] set radius 2.30
    [atomselect $molid "resname DMPA and name C2"] set radius 2.30

    # phosphate atoms in ADP or ATP
    [atomselect $molid "name PA PB PG"] set radius 2.35
    [atomselect $molid "(resname ATP and (name O3A O3B)) or (resname ADP and name O3A)"] set radius 1.65
    [atomselect $molid "(name O1A O2A O1B O2B O1G O2G O3G) or (resname ADP and name O3B)"] set radius 1.49

    # phosphate atoms in phosphotyrosine
    [atomselect $molid "resname TYR and name P1"] set radius 2.35
    [atomselect $molid "resname TYR and (name O2 O3 O4)"] set radius 1.49
}

# **************** For AMBER ****************

# assign radii based on a given AMBER PARM7 file
proc ::cafe::mmpbsa::assign_radii_parm7 { molid args } {
    set rs [parse_parm7 $args]
    [atomselect $molid "all"] set radius $rs
}

# unfortunately, VMD discards the radii info in a PARM7 file,
# so I have to re-parse it myself.
proc ::cafe::mmpbsa::parse_parm7 { fname } {
    if { [catch {set fp [open $fname r]} err] } {
        show -err "Failed to open PARM7 file: $fname"
    }

    set natoms 0
    set nr 0
    set rs { }

    if { [gets $fp line] < 0 } {
        show -err "Found an empty file: $fname"
    } else {
        if { [string range $line 0 7] ne "%VERSION" } {
            show -err "Unknown or broken PARM7 file: $fname"
        }
    }

    while { [gets $fp line] >= 0 } {
        if { [regexp {^%FLAG POINTERS} $line -> s] } {
            # skip the format line
            if { [gets $fp line] < 0 } { break }
            if { [gets $fp line] < 0 } { break }
            set natoms [expr [string range $line 0 7] * 1]
            if { $natoms <= 0 } { break }
        }

        if { [regexp {^%FLAG RADII} $line -> s] } {
            # skip the format line
            if { [gets $fp line] < 0 } { break }
            if { $natoms <= 0 } { break }
            set n [expr $natoms / 5]
            set rest [expr $natoms % 5]

            for { set i 0 } { $i < $n } { incr i } {
                if { [gets $fp line] < 0 } { break }
                for { set j 0 } { $j < 5 } { incr j } {
                    set begin [expr 16 * $j]
                    set end [expr $begin + 15]
                    lappend rs [expr [string range $line $begin $end] * 1.0]
                    incr nr
                }
            }

            if { [gets $fp line] < 0 } { break }

            for { set j 0 } { $j < $rest } { incr j } {
                set begin [expr 16 * $j]
                set end [expr $begin + 15]
                lappend rs [expr [string range $line $begin $end] * 1.0]
                incr nr
            }
            break
        }
    }

    close $fp

    # do not use llength
    if { $natoms <= 0 || $nr != $natoms } {
        show -err "Unknown or broken PARM7 file: $fname"
    }

    return $rs
}

# Swanson et al (J Chem Theory Comput 2007, 3, 170-183)
proc ::cafe::mmpbsa::assign_radii_swanson { molid } {
    show -err "Not implemented yet"
}

# Tan et al (J Phys Chem B 2006, 110, 18680-18687)
proc ::cafe::mmpbsa::assign_radii_tan { molid } {
    show -err "Not implemented yet"
}

# Yamagishi et al (J Comp Chem 2014, 35, 2132-2139)
proc ::cafe::mmpbsa::assign_radii_yamagishi { molid } {
    show -err "Not implemented yet"
}

proc ::cafe::mmpbsa::calc_pb { molid prefix selstr } {
    variable pb
    variable kt2kc
    variable j2cal
    variable first
    variable stride
    variable debug

    if { $pb == 1 } {
        set pb_list [run_delphi $molid $prefix $selstr]
        set conv $kt2kc
    } elseif { $pb == 2 } {
        set pb_list [run_apbs $molid $prefix $selstr]
        set conv $j2cal
    } else {
        show -err "Unknown PB type"
    }

    if { $debug > 0 } {
        set tfmt "#%14s %15s"
        set dfmt "%15d %15.4f"
        set fname ${prefix}_pb.log
        set title [format $tfmt Frame PB]
        set frames { }

        for { set i 0 } { $i < [llength $pb_list] } { incr i } {
            lappend frames [expr $first + $i*$stride]
        }

        set data { }
        lappend data $frames
        lappend data $pb_list

        write_log $fname $title $dfmt $data
    }

    set pb_list [vecscale $conv $pb_list]

    return $pb_list
}

# solve PBE by DelPhi
proc ::cafe::mmpbsa::run_delphi { molid prefix selstr } {
    variable first
    variable stride
    variable pb_exe
    variable pb_siz
    variable pb_crg
    variable pb_indi
    variable pb_exdi
    variable pb_scale
    variable pb_perfil
    variable pb_prbrad
    variable pb_linit
    variable pb_maxc
    variable pb_bndcon
    variable debug

    set sel [atomselect $molid $selstr]

    set pb_list { }

    for { set i 0 } { $i < [molinfo $molid get numframes] } { incr i } {
        molinfo $molid set frame $i
        set iframe [expr $first + $i*$stride]

        set tmp_prefix ${prefix}_pb_tmp_${iframe}

        if { $pb_siz eq "" || $pb_crg eq "" } {
            # DelPhi cannot read files generated by "animate write pqr"
            set tmp_pqr ${tmp_prefix}.pqr
            write_pqr $sel $tmp_pqr
        } else {
            set tmp_pdb ${tmp_prefix}.pdb
            $sel writepdb $tmp_pdb
        }

        foreach { x y z } [measure center $sel] { break }

        set pb_out 0.0
        set tmp_prm ${tmp_prefix}.delphi
        set delphiin [open $tmp_prm w]

        if { $pb_siz eq "" || $pb_crg eq "" } {
            puts $delphiin "in(modpdb4,file=\"$tmp_pqr\",format=pqr)"
        } else {
            puts $delphiin "in(pdb,file=\"$tmp_pdb\")"
            puts $delphiin "in(siz,file=\"$pb_siz\")"
            puts $delphiin "in(crg,file=\"$pb_crg\")"
        }
        puts $delphiin "acenter($x,$y,$z)"
        puts $delphiin "indi=$pb_indi"
        puts $delphiin "exdi=$pb_exdi"
        puts $delphiin "scale=$pb_scale"
        puts $delphiin "perfil=$pb_perfil"
        puts $delphiin "prbrad=$pb_prbrad"
        puts $delphiin "linit=$pb_linit"
        puts $delphiin "maxc=$pb_maxc"
        puts $delphiin "bndcon=$pb_bndcon"
        puts $delphiin "energy(s,g)"

        close $delphiin

        if { [file exists ARCDAT] } { file delete ARCDAT }

        set tmp_log ${tmp_prefix}.log
        exec $pb_exe $tmp_prm > $tmp_log

        # The polar solvation energy is simply C.R.F.Ener(exdi/indi)
        # obtained from the DelPhi log.
        # the unit is kT
        set tmp_out [parse_delphi $tmp_log]
        set pb_out $tmp_out
        lappend pb_list $pb_out

        if { $debug < 2 } { cleanup ${prefix}_pb }
    }

    $sel delete

    if { [file exists ARCDAT] } { file delete ARCDAT }

    return $pb_list
}

# generate a PQR file for DelPhi calcultion
# DelPhi cannot read the PQR written by "animate write pqr"
proc ::cafe::mmpbsa::write_pqr { sel fname } {
    set fp [open $fname w]
    foreach n [$sel get serial] name [$sel get name] \
            rname [$sel get resname] rid [$sel get resid] \
            x [$sel get x] y [$sel get y] z [$sel get z] \
            q [$sel get charge] r [$sel get radius] {
        set name [format "%-3s" $name]
        puts $fp [format "ATOM  %5i %4s %3s  %4i    %8.3f%8.3f%8.3f %7.4f%7.4f" \
                  $n $name $rname $rid $x $y $z $q $r]
    }
    close $fp
}

# generate a ligand size file for DelPhi calcultion
proc ::cafe::mmpbsa::write_siz { molid rname fname } {
    set fp [open $fname w]
    puts $fp atom__res_radius_

    set sel [atomselect $molid "resname $rname"]

    foreach name [$sel get name] r [$sel get radius] {
        puts $fp [format "%-6s%-4s%7.3f" $name $rname $r]
    }

    $sel delete

    close $fp
}

# generate a ligand charge file for DelPhi calcultion
proc ::cafe::mmpbsa::write_crg { molid rname fname } {
    set fp [open $fname w]
    puts $fp atom__resnumbc_charge_

    set sel [atomselect $molid "resname $rname"]

    foreach name [$sel get name] q [$sel get charge] {
        puts $fp [format "%-6s%-8s%8.3f" $name $rname $q]
    }

    $sel delete

    close $fp
}

# parse a DelPhi output file
proc ::cafe::mmpbsa::parse_delphi { fname } {
    if { [catch { set fp [open $fname r] } err] } {
        show -err "Failed to open DelPhi output file: $fname"
    }

    set found 0

    while { [gets $fp line] >= 0 } {
        if { [regexp {^ corrected reaction field energy\s*:(.*)kt} $line -> s] } {
            set s [string trim $s]
            set found 1
            break
        }
    }

    close $fp

    if { !$found } {
        show -err "Unknown or broken DelPhi output file: $fname"
    }

    return $s
}

# solve PBE by APBS
proc ::cafe::mmpbsa::run_apbs { molid prefix selstr } {
    variable first
    variable stride
    variable pb_exe
    variable pb_indi
    variable pb_exdi
    variable pb_scale
    variable pb_perfil
    variable pb_prbrad
    variable pb_bcfl
    variable pb_chgm
    variable pb_srfm
    variable pb_swin
    variable pb_sdens
    variable debug

    set sel [atomselect $molid $selstr]
    set grid [expr 1.0/$pb_scale]
    set pb_list { }

    for { set i 0 } { $i < [molinfo $molid get numframes] } { incr i } {
        molinfo $molid set frame $i
        set iframe [expr $first + $i*$stride]

        set tmp_prefix ${prefix}_pb_tmp_${iframe}
        set tmp_pqr ${tmp_prefix}.pqr
        write_pqr $sel $tmp_pqr

        set x [$sel get x]
        set y [$sel get y]
        set z [$sel get z]
        set r [$sel get radius]

        foreach { xmin xmax } [minmax $x $r] { break }
        foreach { ymin ymax } [minmax $y $r] { break }
        foreach { zmin zmax } [minmax $z $r] { break }

        set len [vecsub [list $xmax $ymax $zmax] [list $xmin $ymin $zmin]]
        set dim { }
        set raw [vecscale [expr 100.0/$pb_perfil*$pb_scale] $len]

        # calculate optimal numbers of grid points for APBS
        # assuming nlev=4, so n=32*c+1
        foreach x $raw {
            set n [expr entier(ceil($x))]
            set r [expr ($n-1)%32]
            if { $r } {
                set c [expr ($n-1)/32+1]
                set n [expr 32*$c+1 ]
            }
            if { $n < 33 } {
                set n 33
            }
            lappend dim $n
        }

        set tmp_prm ${tmp_prefix}.apbs
        set apbsin [open $tmp_prm w]

        puts $apbsin "read"
        puts $apbsin "    mol pqr $tmp_pqr"
        puts $apbsin "end"
        puts $apbsin "elec name pol"
        puts $apbsin "    mg-manual"
        puts $apbsin "    mol 1"
        puts $apbsin "    dime [lindex $dim 0] [lindex $dim 1] [lindex $dim 2]"
        puts $apbsin "    grid $grid $grid $grid"
        puts $apbsin "    gcent mol 1"
        puts $apbsin "    lpbe"
        puts $apbsin "    bcfl $pb_bcfl"
        puts $apbsin "    pdie $pb_indi"
        puts $apbsin "    sdie $pb_exdi"
        puts $apbsin "    chgm $pb_chgm"
        puts $apbsin "    srfm $pb_srfm"
        puts $apbsin "    srad $pb_prbrad"
        puts $apbsin "    swin $pb_swin"
        puts $apbsin "    sdens $pb_sdens"
        puts $apbsin "    temp 298.15"
        puts $apbsin "    calcenergy total"
        puts $apbsin "    calcforce no"
        puts $apbsin "end"
        puts $apbsin "elec name ref"
        puts $apbsin "    mg-manual"
        puts $apbsin "    mol 1"
        puts $apbsin "    dime [lindex $dim 0] [lindex $dim 1] [lindex $dim 2]"
        puts $apbsin "    grid $grid $grid $grid"
        puts $apbsin "    gcent mol 1"
        puts $apbsin "    lpbe"
        puts $apbsin "    bcfl $pb_bcfl"
        puts $apbsin "    pdie $pb_indi"
        puts $apbsin "    sdie $pb_indi"
        puts $apbsin "    chgm $pb_chgm"
        puts $apbsin "    srfm $pb_srfm"
        puts $apbsin "    srad $pb_prbrad"
        puts $apbsin "    swin $pb_swin"
        puts $apbsin "    sdens $pb_sdens"
        puts $apbsin "    temp 298.15"
        puts $apbsin "    calcenergy total"
        puts $apbsin "    calcforce no"
        puts $apbsin "end"
        puts $apbsin "print"
        puts $apbsin "    elecEnergy pol - ref"
        puts $apbsin "end"
        puts $apbsin "quit"

        close $apbsin

        if { [file exists io.mc] } { file delete io.mc }

        set tmp_log ${tmp_prefix}.log
        exec $pb_exe $tmp_prm > $tmp_log

        # the unit is kcal/mol
        set pb_out [parse_apbs $tmp_log -elec]
        lappend pb_list $pb_out

        if { $debug < 2 } { cleanup ${prefix}_pb }
    }

    $sel delete

    if { [file exists io.mc] } { file delete io.mc }

    return $pb_list
}

proc ::cafe::mmpbsa::minmax { x r } {
    set xmin 9999
    set xmax -9999

    foreach n [vecsub $x $r] {
        if { $n < $xmin } { set xmin $n }
    }

    foreach n [vecadd $x $r] {
        if { $n > $xmax } { set xmax $n }
    }

    return [list $xmin $xmax]
}

# parse an APBS output file
proc ::cafe::mmpbsa::parse_apbs { fname opt } {
    if { [catch { set fp [open $fname r] } err] } {
        show -err "Failed to open APBS output file: $fname"
    }

    set found 0

    switch -- $opt {
        -elec {
            while { [gets $fp line] >= 0 } {
                if { [regexp {^  Global net ELEC energy\s*=(.*)kJ/mol} $line -> s] } {
                    set s [string trim $s]
                    set found 1
                    break
                }
            }
        }
        -sasa {
            while { [gets $fp line] >= 0 } {
                if { [regexp {^Total solvent accessible surface area:(.*)A\^2} $line -> s] } {
                    set s [string trim $s]
                    set found 1
                    break
                }
            }
        }
        default { show -err "Typo?" }
    }

    close $fp

    if { !$found } {
        show -err "Unknown or broken APBS output file: $fname"
    }

    return $s
}

# ######################################################################
#                             SA related
# ######################################################################
proc ::cafe::mmpbsa::calc_sa { molid prefix selstr } {
    variable first
    variable stride
    variable sa
    variable sa_prbrad
    variable sa_samples
    variable sa_gamma
    variable sa_beta
    variable debug

    if { $sa == 1 } {
        set sa_list [run_vmd $molid $prefix $selstr]
    } elseif { $sa == 2 } {
        set sa_list [run_apbs_apol $molid $prefix $selstr]
    } else {
        show -err "Unknown non-polar type"
    }

    if { $debug > 0 } {
        set tfmt "#%14s %15s"
        set dfmt "%15d %15.4f"
        set fname ${prefix}_sa.log
        set title [format $tfmt Frame SASA]
        set frames { }

        for { set i 0 } { $i < [llength $sa_list] } { incr i } {
            lappend frames [expr $first + $i*$stride]
        }

        set data { }
        lappend data $frames
        lappend data $sa_list

        write_log $fname $title $dfmt $data
    }

    set offset [lrepeat [llength $sa_list] $sa_beta]
    set sa_list [vecadd [vecscale $sa_gamma $sa_list] $offset]

    return $sa_list
}

# calculate SASA by VMD
proc ::cafe::mmpbsa::run_vmd { molid prefix selstr } {
    variable stride
    variable sa_prbrad
    variable sa_samples
    variable sa_gamma
    variable sa_beta

    set sa_list { }
    set sel [atomselect $molid $selstr]

    for { set i 0 } { $i < [molinfo $molid get numframes] } { incr i } {
        molinfo $molid set frame $i

        # the unit is angstrom^2
        set val [measure sasa $sa_prbrad $sel -samples $sa_samples]
        lappend sa_list $val
    }

    $sel delete

    return $sa_list
}

# calculate SASA/SAV/WCA by APBS
proc ::cafe::mmpbsa::run_apbs_apol { molid prefix selstr } {
    variable first
    variable stride
    variable sa_exe
    variable sa_prbrad
    variable sa_samples
    variable sa_gamma
    variable sa_beta
    variable debug

    set sa_list { }
    set sel [atomselect $molid $selstr]

    for { set i 0 } { $i < [molinfo $molid get numframes] } { incr i } {
        molinfo $molid set frame $i
        set iframe [expr $first + $i*$stride]

        set tmp_prefix ${prefix}_sa_tmp_${iframe}
        set tmp_pqr ${tmp_prefix}.pqr
        write_pqr $sel $tmp_pqr

        set tmp_prm ${tmp_prefix}.apbs
        set apbsin [open $tmp_prm w]

        puts $apbsin "read"
        puts $apbsin "    mol pqr $tmp_pqr"
        puts $apbsin "end"
        puts $apbsin "apolar name apol"
        puts $apbsin "    mol 1"
        puts $apbsin "    gamma 1"
        puts $apbsin "    press 0"
        puts $apbsin "    bconc 0"
        puts $apbsin "    dpos 0.01"
        puts $apbsin "    grid 0.5 0.5 0.5"
        puts $apbsin "    sdens $sa_samples"
        puts $apbsin "    srad $sa_prbrad"
        puts $apbsin "    srfm sacc"
        puts $apbsin "    swin 0.3"
        puts $apbsin "    temp 298.15"
        puts $apbsin "    calcenergy total"
        puts $apbsin "    calcforce no"
        puts $apbsin "end"
        puts $apbsin "print"
        puts $apbsin "    apolEnergy apol"
        puts $apbsin "end"
        puts $apbsin "quit"

        close $apbsin

        if { [file exists io.mc] } { file delete io.mc }

        set tmp_log ${tmp_prefix}.log
        exec $sa_exe $tmp_prm > $tmp_log

        # the unit is angstrom^2
        set sa_out [parse_apbs $tmp_log -sasa]
        lappend sa_list $sa_out

        if { $debug < 2 } { cleanup ${prefix}_sa }
    }

    $sel delete

    if { [file exists io.mc] } { file delete io.mc }

    return $sa_list
}

