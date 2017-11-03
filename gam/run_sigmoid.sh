#!/bin/bash
/home/matlab/matlab_bin/bin/matlab -nodisplay -nodesktop -nosplash -r \
      "try, run('run_sigmoid.m'), catch me, fprintf('%s / %s\n',me.identifier,me.message), exit(1), end, exit(0);"
      echo "matlab exit code: $?"
