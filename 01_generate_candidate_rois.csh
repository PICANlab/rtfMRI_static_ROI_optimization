#!/usr/bin/csh

# Generate 810 candidate 7-mm spherical ROIs for each neurofeedback target.

set cimg = $1
set R = 7

mkdir -p amygballs
mkdir -p hipsballs


# ----------------------------------------------------------------------
# Amygdala
#
# 9 x 10 x 9 = 810 candidate ROIs
# Original AFNI-space center: 21 5 -16
# ----------------------------------------------------------------------

foreach px ('17' '18' '19' '20' '21' '22' '23' '24' '25')
    foreach py ('1' '2' '3' '4' '5' '6' '7' '8' '9' '10')
        foreach pz ('-12' '-13' '-14' '-15' '-16' '-17' '-18' '-19' '-20')

            3dcalc \
                -a $cimg \
                -prefix amygballs/amygball_${px}_${py}_${pz} \
                -expr "1*step($R*$R-(x-$px)*(x-$px)-(y-$py)*(y-$py)-(z-$pz)*(z-$pz))"

        end
    end
end


# ----------------------------------------------------------------------
# Parietal (HIPS)
#
# 9 x 10 x 9 = 810 candidate ROIs
# Original AFNI-space center: 42 48 48
# ----------------------------------------------------------------------

foreach px ('38' '39' '40' '41' '42' '43' '44' '45' '46')
    foreach py ('44' '45' '46' '47' '48' '49' '50' '51' '52' '53')
        foreach pz ('44' '45' '46' '47' '48' '49' '50' '51' '52')

            3dcalc \
                -a $cimg \
                -prefix hipsballs/hipsball_${px}_${py}_${pz} \
                -expr "1*step($R*$R-(x-$px)*(x-$px)-(y-$py)*(y-$py)-(z-$pz)*(z-$pz))"

        end
    end
end
