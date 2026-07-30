# amr_classify_scaled.py
# Pure Python Classifier with Feature Scaling (Min-Max Normalization)

import csv
import math
from collections import defaultdict

print("[1/5] Loading generated feature dataset...")

data = []
with open('amr_protein_features.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append(row)

ignore_cols = {'ID', 'Family', 'Sequence'}
feature_cols = [col for col in data[0].keys() if col not in ignore_cols]

print(f"Loaded {len(data)} samples with {len(feature_cols)} features each.")

print("[2/5] Applying Min-Max Feature Scaling (Normalizing to [0, 1])...")

# Find Min and Max for each feature
min_vals = {col: float('inf') for col in feature_cols}
max_vals = {col: float('-inf') for col in feature_cols}

for sample in data:
    for col in feature_cols:
        val = float(sample[col])
        if val < min_vals[col]:
            min_vals[col] = val
        if val > max_vals[col]:
            max_vals[col] = val

# Normalize dataset: (value - min) / (max - min)
scaled_data = []
for sample in data:
    scaled_sample = {'ID': sample['ID'], 'Family': sample['Family']}
    for col in feature_cols:
        val = float(sample[col])
        denom = (max_vals[col] - min_vals[col])
        scaled_sample[col] = (val - min_vals[col]) / denom if denom > 0 else 0.0
    scaled_data.append(scaled_sample)

print("[3/5] Splitting data into Train and Test sets...")

# Stratified-like stride shuffle
data_shuffled = []
for i in range(len(scaled_data)):
    data_shuffled.append(scaled_data[(i * 37) % len(scaled_data)])

split_idx = int(len(data_shuffled) * 0.8)
train_set = data_shuffled[:split_idx]
test_set = data_shuffled[split_idx:]

print(f"Train size: {len(train_set)} | Test size: {len(test_set)}")

print("[4/5] Computing Feature Centroids for each Gene Family...")

centroids = defaultdict(lambda: [0.0] * len(feature_cols))
counts = defaultdict(int)

for sample in train_set:
    fam = sample['Family']
    counts[fam] += 1
    for i, feat in enumerate(feature_cols):
        centroids[fam][i] += sample[feat]

for fam in centroids:
    for i in range(len(feature_cols)):
        centroids[fam][i] /= counts[fam]

def euclidean_distance(vec1, vec2):
    return math.sqrt(sum((v1 - v2) ** 2 for v1, v2 in zip(vec1, vec2)))

print("[5/5] Evaluating Scaled Model Accuracy...\n")

correct = 0
family_correct = defaultdict(int)
family_total = defaultdict(int)

for sample in test_set:
    true_fam = sample['Family']
    sample_vec = [sample[feat] for feat in feature_cols]
    
    best_fam = None
    min_dist = float('inf')
    
    for fam, centroid in centroids.items():
        dist = euclidean_distance(sample_vec, centroid)
        if dist < min_dist:
            min_dist = dist
            best_fam = fam
            
    family_total[true_fam] += 1
    if best_fam == true_fam:
        correct += 1
        family_correct[true_fam] += 1

total_accuracy = (correct / len(test_set)) * 100

print("=" * 55)
print(f"{'SCALED PREDICTION EVALUATION REPORT':^55}")
print("=" * 55)
print(f"{'FAMILY':<10} | {'CORRECT':<10} | {'TOTAL':<10} | {'ACCURACY':<10}")
print("-" * 55)
for fam in sorted(counts.keys()):
    c = family_correct[fam]
    t = family_total[fam]
    acc = (c / t * 100) if t > 0 else 0
    print(f"{fam:<10} | {c:<10} | {t:<10} | {acc:.1f}%")
print("=" * 55)
print(f"Overall Scaled Test Accuracy: {total_accuracy:.2f}%\n")