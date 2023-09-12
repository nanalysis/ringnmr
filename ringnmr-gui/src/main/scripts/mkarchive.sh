#!/bin/zsh

JREHOME='/Users/brucejohnson/Development/jres'
jversion='jdk-11.0.9.1+1-jre'
jversion='jdk-16.0.51-jre'
jversion='zulu17.40.19-ca-fx-jre17.0.6'

#zulu17.44.15-ca-fx-jre17.0.8-win_x64
#zulu17.44.15-ca-fx-jre17.0.8-linux_x64
#zulu17.44.15_1-ca-fx-jre17.0.8-macosx_aarch64

jversionMac='zulu17.44.15_1-ca-fx-jre17.0.8'
jversion='zulu17.44.15-ca-fx-jre17.0.8'

PRGSCRIPT=ringnmr-gui

dir=`pwd`
PRG="$(basename $dir)"

if [ -e "installers" ]
then
    rm -rf installers
fi

#for os in "linux-amd64"
#for os in "macosx_x64"
for os in "macosx_aarch64" "macosx_x64" "linux_x64" "win_x64"
do
    if [[ $os == "linux_x64" ]]
    then
        jreFileName=${jversion}-${os}
    elif [[ $os == "macosx_x64" ]]
    then
        jreFileName=${jversionMac}-${os}
    elif [[ $os == "macosx_aarch64" ]]
    then
        jreFileName=${jversionMac}-${os}
    else
        jreFileName=${jversion}-${os}
    fi
    echo $jreFileName

    dir=installers/$os
    if [ -e $dir ]
    then
         rm -rf $dir
    fi

    mkdir -p $dir
    cd $dir
    cp -r -p ../../target/${PRG}-*-bin/${PRG}* .
    sdir=`ls -d ${PRG}-*`
    cd $sdir
    echo $sdir
    cp -r -p ../../../../ringnmr/target/*-bin/ring*/ringnmr .
    cp -r -p ../../../../ringnmr/target/*-bin/ring*/ringnmr.bat .

    rm lib/javafx*

    if [[ $os == "linux_x64" ]]
    then
        rm lib/*-mac*
        rm lib/*-win*
    fi

    if [[ $os == "win_x64" ]]
    then
        rm lib/*-linux*
        rm lib/*-mac*
    fi

    if [[ $os == "macosx_x64" ]]
    then
        cp -R -p ${JREHOME}/$jreFileName .
        rm lib/*-linux*
        rm lib/*-win*
        rm lib/*-macosx_aa**
    elif [[ $os == "macosx_aarch64" ]]
    then
        cp -R -p ${JREHOME}/$jreFileName .
        rm lib/*-linux*
        rm lib/*-win*
        rm lib/*-macosx_x8**
    else
        cp -r -p ${JREHOME}/$jreFileName jre
    fi
    cd ..

    fname=`echo $sdir | tr '.' '_'`
    if [[ $os == "linux_x64" ]]
    then
        tar czf ${fname}_${os}.tar.gz $sdir
    elif [[ $os == "macosx_x64" ]]
    then
        #xattr -cr $sdir
        tar czf ${fname}_${os}.tar.gz $sdir
    else
        zip -r ${fname}_${os}.zip $sdir
    fi
    cd ../..
done
