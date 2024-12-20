mkdir -p ~/project/result/model_fitting/
mkdir -p ~/project/result/model_final/

#### Run test models ####
## 2 logs as ran LOO models on second instance
nohup Rscript ~/project/scripts/models/model_fit_test.R >> ~/project/result/logs/model_fit_test.log 2>&1 &

#### Run GMMAT full models ####
## 3 model_full logs b/c ran out of RAM (1) and had a typo (2)
nohup Rscript ~/project/scripts/models/model_full.R >> ~/project/result/logs/model_full.log 2>&1 &
nohup Rscript ~/project/scripts/models/model_comp.R >> ~/project/result/logs/model_comp.log 2>&1 &

#### Run GMMAT cov check models ####
nohup Rscript ~/project/scripts/models/model_signif_cov.R >> ~/project/result/logs/model_signif_cov.log 2>&1 &
nohup Rscript ~/project/scripts/models/model_loo.R >> ~/project/result/logs/model_loo.log 2>&1 &

aws s3 sync ~/project/result s3://hawn-reporttb-results2
# aws s3 sync s3://hawn-reporttb-results2 ~/project/result 
