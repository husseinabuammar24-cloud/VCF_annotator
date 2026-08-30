# VCF_annotator



Note :
The error above (grayed out):

python3 variant_annotator
...can't find '__main__' module in '.../VCF_annotator'

That happened because python3 <folder> only works if the folder has an __init__.py/__main__.py making it a runnable package. variant_annotator doesn't have one, so Python couldn't figure out what to execute. Not a bug in your code — just the wrong invocation.



# Engineering & Bioinformatics Lessons Learned

### 1. Versioning and Reference Context
* Genomic coordinates are version-dependent.
* Structural shifts happen between genome updates.
* Missing assembly metadata breaks validation pipelines.
* Always check the underlying reference schema first.
* Embed version tracking directly into data flows.

### 2. Interface Layer Selection
* Application sites are not API systems.
* Web domains return unparseable HTML text.
* Target dedicated developer subdomains instead.
* Choose `POST` methods for complex arrays.
* Offload orientation and strand logic to endpoints.

### 3. Classification of Error Types
* Separate flaky errors from design bugs.
* Identical inputs producing identical errors mean faulty logic.
* Changing answers across identical runs mean flaky infrastructure.
* Deploy automated retry loops for transient infrastructure.
* Refuse to retry structural code errors.

### 4. Diagnostic Logging Fidelity
* Catching exceptions with generic fallbacks hides bugs.
* Do not sacrifice data clarity for console aesthetics.
* Propagate literal HTTP response codes immediately.
* Print raw backend server error string payloads.
* Trace true errors to skip hours of guessing.

### 5. Mock Data Integrity
* Hand-written test files cause artificial data bugs.
* Human blind spots create invalid biological mockups.
* Ground test suites in validated reference data.
* Match experimental profiles against verified ground truth.

### 6. Code Scope and Scaling
* Misplaced indentation layers break variable scopes.
* Trace active variable lifetimes across execution blocks.
* Space connection queries to prevent network blocks.
* Transition from singular loops to batch arrays early.
Use code with caution.




                ┌──────────────────────┐
                │      VCF file        │
                │ real_grch38_test.vcf │
                └──────────┬───────────┘
                           │
                           ▼
                ┌──────────────────────┐
                │      parse_vcf()     │
                │                      │
                │ VCF text → Variant   │
                │ objects              │
                └──────────┬───────────┘
                           │
                           ▼
                ┌──────────────────────┐
                │ split_multiallelic() │
                │                      │
                │ One ALT per Variant  │
                └──────────┬───────────┘
                           │
                           ▼
              ┌───────────────────────────┐
              │ determine_variant_impact()│
              │                           │
              │ 1. Clean chr name         │
              │ 2. Handle indel notation  │
              │ 3. Build VEP request      │
              │ 4. Send HTTP request      │
              │ 5. Retry 5xx errors       │
              └────────────┬──────────────┘
                           │
                           ▼
                 ┌─────────────────────┐
                 │   Ensembl VEP REST  │
                 │      GRCh38         │
                 └─────────┬───────────┘
                           │
                           ▼
                 ┌─────────────────────┐
                 │   VEP JSON result   │
                 │                     │
                 │ transcript_consequences
                 │ consequence_terms   │
                 └─────────┬───────────┘
                           │
                           ▼
                ┌──────────────────────┐
                │ Impact classification│
                │                      │
                │ frameshift →         │
                │ Frameshift/Indel     │
                │ stop_gained →        │
                │ Nonsense             │
                │ missense → Missense  │
                │ synonymous →         │
                │ Synonymous           │
                └──────────┬───────────┘
                           │
                           ▼
                ┌──────────────────────┐
                │      Terminal        │
                │                      │
                │ chr19:... T>C | SNV  │
                │              |Missense│
                └──────────────────────┘