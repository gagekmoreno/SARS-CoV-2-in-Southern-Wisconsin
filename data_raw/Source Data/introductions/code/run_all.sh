python3 code/infer_tree.py
python3 code/infer_clock_mugration.py
for i in $(seq 1 100); do; python3 code/estimate_importations.py --bootstrap_replicate $i; done; 
python3 code/rarefaction.py --n_samples 20 40 60 80 100 120 140 160 180 200 220 240
python3 code/analyze_results.py