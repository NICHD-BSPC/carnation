#### Load New or Edit Existing Data
-----------------------------------

Interface for importing and processing new RNA-Seq datasets into Carnation or to edit
an existing Carnation dataset.

**What it does:**
- Guides you through the process of importing count data and metadata
- Processes differential expression results
- Creates a Carnation object that can be saved and reused

**Required files:**
- **Counts file**: Tab-delimited text file with gene IDs in the first column and sample counts in subsequent columns
- **Metadata file**: Tab-delimited text file with sample names in the first column and condition information in subsequent columns
- **DE results file(s)**: Tab-delimited text files from differential expression analysis tools (currently supports *DESeq2* format)

**Optional files:**
- **Functional enrichment results**: Tab-delimited text files from enrichment analysis tools (currently supports *clusterProfiler* output)

**How to use:**
1. Click "Add counts" to upload count matrix and metadata files
2. Click "Add DE results" to upload differential expression results
3. Optionally add functional enrichment results
4. Provide an analysis label and select a save location
5. Click "Create object" to process and save the data

**Important settings:**
- **Compress RDS**: Toggle compression of the saved file (uncompressed files load faster but are larger)
- **Force overwrite**: Enable to replace existing files with the same name
- **User group**: Select which user group can access this dataset

**Note:** Large datasets may take several minutes to process. The application will reload automatically once the data is successfully imported.
