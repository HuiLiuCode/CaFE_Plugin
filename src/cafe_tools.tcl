package provide cafe_tools 1.0

namespace eval ::cafe::tools:: {
    namespace export show timer check_file check_exe check_string check_bool \
                     check_int check_pos_int check_nneg_int check_real \
                     check_pos_real check_nneg_real write_item write_log \
                     parse_namd cleanup calc_stats
}

# print info/error messages
proc ::cafe::tools::show { opt { msg "" } } {
    switch -- $opt {
        -info {
            puts "CaFE) $msg"
        }
        -err {
            error "CaFE) Error: $msg"
        }
        default {
            error ""
        }
    }
}

# a timer
proc ::cafe::tools::timer { start } {
    set tmp [expr [clock seconds] - $start]
    set day [expr $tmp / 86400]
    set tmp [expr $tmp % 86400]
    set hrs [expr $tmp / 3600]
    set tmp [expr $tmp % 3600]
    set min [expr $tmp / 60]
    set sec [expr $tmp % 60]
    return [list $day $hrs $min $sec]
}

proc ::cafe::tools::check_file { s n } {
    if { [file exists $s] } {
        return $s
    } else {
        show -err "Cannot find file '$s' for '$n'"
    }
}

proc ::cafe::tools::check_exe { s n } {
    if { [file executable $s] } {
        return $s
    } else {
        # $tcl_platform(pathSeparator) is only in Tcl 8.6 and later
        # do nothing for macintosh and others
        switch -- $::tcl_platform(platform) {
            -unix { set sep ":" }
            -windows { set sep ";" }
            default { return $s }
        }

        foreach path [split $::env(PATH) $sep] {
            set exe [file join $path $s]
            if [file executable $exe] { return $s }
        }
        show -err "Expected an executable file for '$n' but got '$s'"
    }
}

proc ::cafe::tools::check_string { s n } {
    if { $s ne "" } {
        return [string tolower $s]
    } else {
        show -err "Expected a non-empty string for '$n'"
    }
}

proc ::cafe::tools::check_bool { s n } {
    if { [string is boolean $s] } {
        return $s
    } else {
        show -err "Expected a boolean for '$n' but got '$s'"
    }
}

proc ::cafe::tools::check_int { s n } {
    if { [string is integer -strict $s] } {
        return $s
    } else {
        show -err "Expected an integer for '$n' but got '$s'"
    }
}

proc ::cafe::tools::check_pos_int { s n } {
    if { [string is integer -strict $s] && $s > 0 } {
        return $s
    } else {
        show -err "Expected a positive integer for '$n' but got '$s'"
    }
}

proc ::cafe::tools::check_nneg_int { s n } {
    if { [string is integer -strict $s] && $s >= 0 } {
        return $s
    } else {
        show -err "Expected a nonnegtive integer for '$n' but got '$s'"
    }
}

proc ::cafe::tools::check_real { s n } {
    if { [string is double -strict $s] } {
        return $s
    } else {
        show -err "Expected a real for '$n' but got '$s'"
    }
}

proc ::cafe::tools::check_pos_real { s n } {
    if { [string is double -strict $s] && $s > 0 } {
        return $s
    } else {
        show -err "Expected a positive real for '$n' but got '$s'"
    }
}

proc ::cafe::tools::check_nneg_real { s n } {
    if { [string is double -strict $s] && $s >= 0 } {
        return $s
    } else {
        show -err "Expected a non-negtive real for '$n' but got '$s'"
    }
}

proc ::cafe::tools::write_item { fp name data } {
    set fmt "   %-6s %15d %15.4f %15.4f"
    foreach { mean sd } [calc_stats $data] { break }
    set n [llength $data]
    puts $fp [format $fmt $name $n $mean $sd]
}

proc ::cafe::tools::write_log { fname title fmt data } {
    if { [catch { set log [open $fname w] } err] } {
        show -err "Failed to write log file: $fname"
    }

    puts $log $title

    for { set i 0 } { $i < [llength [lindex $data 0]] } { incr i } {
        set s "\"$fmt\""
        for { set j 0 } { $j < [llength $data] } { incr j } {
            append s " [lindex [lindex $data $j] $i]"
        }
        puts $log [eval format $s]
    }

    close $log
}

# parse a NAMD output file
proc ::cafe::tools::parse_namd { fname } {
    if { [catch { set fp [open $fname r] } err] } {
        show -err "Failed to open NAMD log file: $fname"
    }

    set enerlists { }

    while { [gets $fp line] >= 0 } {
        if { [regexp {^ENERGY:} $line] } {
            set enerlist { }
            foreach item [split $line] {
                if { $item ne "" } {
                    lappend enerlist $item
                }
            }
            # bond angle dihedral improper elec vdw
            lappend enerlists [lrange $enerlist 2 7]
        }
    }

    close $fp

    return $enerlists
}

# remove tmp files
proc ::cafe::tools::cleanup { prefix } {
    foreach f [glob ${prefix}_tmp*] {
        file delete $f
    }
}

# statistical analysis
proc ::cafe::tools::calc_stats { x } {
    # mean:
    #   av = sum(x)/n
    # standard diviation (SD):
    #   sd = sqrt(sum(x-av)^2/n)
    # standard error of the mean (SEM):
    #   sem = sd/sqrt(n)
    set av [vecmean $x]
    set sd [vecstddev $x]

    return [list $av $sd]
}

