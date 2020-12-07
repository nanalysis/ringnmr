@echo off

rem nvjp [ script  [ arg ... ] ]
rem 
rem optional environment variables:
rem
rem JAVA_HOME  - directory of JDK/JRE, if not set then 'java' must be found on PATH
rem CLASSPATH  - colon separated list of additional jar files & class directories
rem JAVA_OPTS  - list of JVM options, e.g. "-Xmx256m -Dfoo=bar"
rem TCLLIBPATH - space separated list of Tcl library directories
rem


if "%OS%" == "Windows_NT" setlocal

set comdnmrver=${project.version}
set comdnmrmain=org.comdnmr.gui.NMRApp


set dir=%~dp0
set javaexe=java

set cp="%dir%ringnmr-gui-%nvjver%.jar;%dir%lib/Manifest.jar"

set testjava=%dir%jre\bin\java.exe

if exist %testjava% (
    set javaexe="%testjava%"
    set cp="%dir%lib/ringnmr-gui-%comdnmrver%.jar;%dir%lib/%Manifest.jar"
)


%javaexe%  -cp %cp% %JAVA_OPTS% %comdnmrmain% %*

