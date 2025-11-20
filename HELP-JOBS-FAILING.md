# HELP! My jobs keep failing!
1. Check your paths and variables
2. Make sure any input .txt, .tsv, .bed, or other files were saved with Unix line endings
3. Check your slurm output file for any errors ```less slurm-<jobID>.out```. These are located in the folder that you ran your code in. For instance, I always run my code out of my "scripts" directory so that I can always find the corresponding log file for the job. The SLURM scheduler uses the following nomenclature: slurm-###.out (e.g. job 5063360 should automatically have slurm-5063360.out.)
   - HELP! There was no slurm-###.out file written: 1) Double check you are looking in the right directory, 2) If you are sure a .out was not written, check to see that the job is/is not running (```squeue --me```), 3) If the job is running (often indefinitely without finishing) and there is no .out, it is possible that your job was assigned to a failing or draining node and the job will need to be restarted. I recommend emailing the server managers with your problem and the job ID to see if they have anything on their end.
4. See if you can understand the error in the slurm-###.out file and address the problem. If you are still stuck, you can always ask ChatGPT (or other bot of choice) to explain what the error code was and how to fix it.
5. If this is your nth error of the day, you have been staring at code all day/week, you can't figure it out even with these steps, and you don't have an immediate deadline, I suggest taking a break, getting a good night's rest, and attacking the error with fresh eyes on a new day. I have resolved many an error after I have left it and come back to it the next day.
6. When in doubt, the answer is 42.
     
