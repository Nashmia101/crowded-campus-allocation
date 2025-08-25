# Crowded Campus Allocation & Typo Detection
Implements a flow network with lower/upper bounds and Ford-Fulkerson to allocate students into classes under time-slot, capacity, and satisfaction constraints, plus a Trie-based word-checking algorithm to detect substitution-only Levenshtein distance = 1 cases efficiently
---

## Problem 1 – A Crowded Campus
- **Task:** Allocate `n` students into `m` proposed classes.  
- **Constraints:**  
  - Each student must be allocated to exactly one class.  
  - Each class must meet its **minimum and maximum occupancy**.  
  - At least `minimumSatisfaction` students must be allocated to one of their **top-5 preferred time slots**.  

### Approach
- Modeled as a **flow network with demands and lower bounds**.  
- Built a graph with **students → timeslots → classes → sink**.  
- Eliminated lower bounds and adjusted demands with a **super source/sink construction**.  
- Applied **Ford-Fulkerson (with BFS augmenting paths)** to check feasibility and compute allocations.  
- Included a **pre-processing phase** to pre-assign students into top-5 slots for satisfaction guarantees.  

---

## Problem 2 – Typo
- **Task:** Identify words from a dictionary that are **exactly one substitution away** from a given suspicious word.  
- **Constraints:** Only substitutions are allowed (no insertions/deletions).  
- **Input:** A list of unique dictionary words and a suspicious word.  
- **Output:** All dictionary words with Levenshtein distance = 1.  

### Approach
- Built a **Trie-based data structure** to efficiently store all words.  
- Implemented a recursive character-by-character traversal that checks substitution differences.  
- Ensured efficient complexity:  
  - **Build:** O(C) where C = total characters.  
  - **Check:** O(J × N) in worst-case, with early stopping on >1 mismatch.  

---

## Results
- Implemented a working **student allocation system** that respects occupancy + satisfaction constraints.  
- Implemented a **string similarity checker** that returns near-miss words efficiently.  
- Both solutions documented with **time & space complexity analysis**.  
  

---

## Author
**Nashmia Shakeel** 
