#!/usr/bin/perl

open(fp,">model/Aux/profasi_version.cc");
print fp "#include \"profasi_version.hh\"\n";
print fp "std::string profasi_version() { return \"1.4.9\"; }\n";
print fp "std::string profasi_git_version() { return \"git:";
print fp `git log -1 | head -n 1 | cut -f2 -d ' ' | tr -d '\n'`;
print fp "\"; }\n";
close(fp);

