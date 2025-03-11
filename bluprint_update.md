<!-- # Implementation Plan

## Documentation, Error and Logging, Visualization,Testing, CI/CD

### 1. Documentation and Tutorials
- **User Guide**: Provide a comprehensive user manual explaining how to install, configure, and run the pipeline, especially in centralized and federated settings.
- **API Documentation**: Document core functions and classes in detail, particularly those interacting with the Flower framework, explaining input/output formats, parameters, and example use cases.
- **Tutorials**: Provide step-by-step tutorials for different use cases, covering how to use the pipeline in federated and centralized settings, including how to simulate distributed environments in the Flower framework.

### 2. Error Handling and Logging
- **Robust Error Handling**: Handle edge cases, missing data, incompatible inputs (such as different SNP formats), and distributed errors that may occur in federated environments (such as client communication failures).
- **Logging System**: Implement logging to track execution processes, including data processing steps, errors, model updates, and performance metrics, especially in federated settings to track data flow and updates from different clients.

### 3. Visualization Tools and Collaboration Features
- **Data Visualization**: Develop visualization tools for key analysis steps, such as quality control metrics, association analysis results, and kinship coefficient matrices. Visualize the global model and client contributions in federated learning in a user-friendly way.
- **Federated Learning Visualization**: Create a visualization interface for the federated learning process, showing each client's contribution to the global model, tracking performance in each iteration, and supporting viewing of model update processes and synchronization status.
- **Collaboration Dashboard**: Build a web-based collaboration dashboard to display federated cooperation progress, showing each client's contributions, model progress, and analysis status.
- **Version Control**: Implement version control for models and analysis scripts to ensure reproducibility between different teams in federated settings.

### 4. Automated Testing and Continuous Integration
- **Unit Testing**: Develop unit tests to ensure the correctness of each module and function in the pipeline, covering both centralized and federated implementations.
- **Continuous Integration (CI)**: Set up continuous integration pipelines to automatically test code changes, ensuring quality control before each release, especially for distributed training and federated learning validation in distributed environments.

## Benchmark and Validation

### 1. Compatibility with External Tools
- **External Tool Integration**: Ensure the pipeline can seamlessly integrate with existing tools (such as PLINK, KING), especially supporting interaction with the Flower framework (relying on PLINK for basic computations).
- **Dataset Integration**: Ensure the pipeline can seamlessly integrate with existing datasets like AllofUs and 1000 Genome Project.

### 2. Model Evaluation and Validation
- **Cross-Validation**: Implement cross-validation techniques to evaluate the robustness of association analysis and kinship inference results, especially in distributed environments.
- **Biological Validation**: Provide integration with tools or databases for functional annotation of SNPs, evaluating the biological relevance of identified associations.
- **Comparative Benchmarking**: Compare performance and results with other common GWAS tools (such as PLINK) on real datasets (federation not required).

### 3. Performance Benchmarking and Testing
- **Performance Benchmarks**: Implement performance benchmarking tools to measure pipeline runtime, memory usage, and scalability in federated settings. Special focus on performance under the Flower framework (on simulated datasets).
- **Comparative Benchmarking**: Compare performance and results with other common GWAS tools (such as PLINK, GEMMA) to ensure excellence in both accuracy and efficiency in federated environments.

## Others

### Security and Privacy Enhancements

- **Channel Encryption**: Implement encryption functionality for data in transit, ensuring security in federated learning and preventing security threats such as man-in-the-middle attacks (check Flower documentation).
- **Access Control**: Establish appropriate user authentication and authorization protocols to manage access to the federated system, ensuring only authorized researchers can access and contribute data.

### Data Transfer Optimization
- **Large-Scale Data Transfer**: Implement optimized data transfer strategies, including compression algorithms to reduce data transfer volume, minimize communication overhead in federated environments, and ensure effective data transmission (check Flower documentation).
- **Data Storage and Synchronization**: Optimize data storage and synchronization strategies to ensure efficiency and synchronicity of data transfer between clients, reducing latency.

### Parallel Computing
- **Parallel Computation**: Ensure support for parallel computing, especially in federated environments, with reasonable scheduling of computational tasks to improve computational efficiency. -->

# Project Timeline (3.10 -  4.30)

## Phase 1: Core Functionality Optimization and Extension - 4 weeks

 - **Implement and Refine Centralized Pipeline Functionality - Assign: @sonam @sitao (1 week)**
  - Finalize API (parameters, return) for each functions in data cleaning (missing values, coin data detection, Hardy-Weinberg equilibrium tests, etc.), association analysis (linear regression, logistic regression), kinship analysis (KING coefficient calculation). and data validation modules to ensure complete functionality in centralized environments.
  - Modulize Quality Control (QC), association analysis, kinship inference, and data validation modules => encapsulate functions in Class (keep state of internal results and outputs in class)
  - Integrate PLINK for basic computations and ensure compatibility with external tools into functions.
  - Implement test functions for each function and module in the pipeline and test correctness and consistency.

- **Dataset Integration and API (I/O) Implementation - Assign: @sonam @xinyue (0.5-1 week)**
  - Implement function/API for loading Human 1000 Genome datasets and test the correctness of the data loading process.
  - Test GWAS functions after reading Human 1000 Genome datasets.
  - Implement function/API for loading AllofUs datasets and test the correctness of the data loading process in AllofUs platform.

- **Implement and Refine Federated GWAS module and pipeline - Assign: @sitao and @sonam (1-2 week)**
  - Refine federated GWAS API
    - Use the **Flower** framework to configure a distributed learning environment and ensure models can run on multiple clients with synchronization and updates.
    - Finalize API parameters and return values
    - Ensure API calls centralized functions
    - Implement test function for each API and test the correctness and consistency
    - Complete simple federated learning training tasks to ensure basic client communication and model synchronization work.
  - Implement single-node running API (command line interface + parameters) for fedGWAS and test the pipeline in a single-node environment.
  - Implement multiple-node running API (commnad line interface + parameters) for fedGWAS and test the pipeline in a multiple-node environment.

  - **Optimize Error Handling, Output Report format and Logging System - Assign: @xinyue @sonam (1 week)**
  - Add comprehensive error handling mechanisms to the pipeline to handle client communication errors, data inconsistencies, etc. in distributed environments.
  - Design which information to save, output, logging, and using **logging** module to create a logging system to record model updates, data synchronization, runtime, etc. for each client.
  - Implement results saving and report generation functionality to store and display quality control metrics, association analysis results, and kinship coefficient matrices.
  - Implement a logging system to track execution processes, including data processing steps, errors, model updates, and performance metrics, especially in federated settings to track data flow and updates from different clients.


### Phase 2: Model Evaluation and Benchmarking - 2 weeks

**Objective:** Enhance model evaluation functionality to ensure result reliability.

- **Implement Cross-Validation Techniques and Biological Validation Processes- Assign: @sonam @xinyue (1 week)**
  - Conduct result validation on each client.
  - Perform biological validation in external tools to ensure SNPs are associated with known diseases or phenotypes.

- **Conduct Performance Benchmarking - Assign: @xinyue @sonam @sitao (1 week)**
  - Design benchmarking plan and experiments details.
  - Run performance benchmarks on large-scale datasets to ensure the pipeline runs efficiently in federated environments.
  - Compare execution time and memory usage under different data scales, especially efficiency in federated environments.
  - **Tech Stack**: Python, cProfile, timeit

### Phase 4: Visualization and Collaboration Features - 1 week
**Objective:** Develop visualization tools and collaboration platforms to enhance user experience and collaboration efficiency.

- **Develop Data Visualization Tools using Tensorboard - @sitao @sonam (1 week)**
  - Develop visualization tools for displaying quality control metrics, association analysis results, and kinship coefficient matrices.

- **Build Web Collaboration Dashboard - TBD**
  - Use **Flask** or **Django** to build a web dashboard displaying the progress, contributions, and results of each client in the federated learning process.

### Phase 3: Tool Documentation, Testing and Report - 1 week

- **Complete User Guide and API Documentation - Assign: @sitao @xinyue**
  - Ensure detailed API documentation, including input/output formats, functionality descriptions, and examples.
  - Write simple tutorials explaining how to set up and use the pipeline.
  - Github page for the project, including installation instructions, usage examples, and API documentation.
  - Report and Paper writing for the project.

<!-- ### Phase 3: Security and Privacy Protection
**Objective:** Enhance data protection and privacy to ensure compliance with relevant security standards.

- **Implement Data Encryption and Access Control Mechanisms:**
  - Data encryption and decryption functionality to ensure data is encrypted during transmission.
  - Implement data compression algorithms to reduce time and bandwidth consumption for data transfer between clients and servers.
  - Use **OAuth2** or **JWT** to implement access control, ensuring only authorized users can access the system.

- **Complete Security Testing for Federated Environments:**
  - Conduct security testing on federated learning environments to ensure client data is not leaked. -->

<!-- ### Phase 5: Cloud Platform and Parallel Computing Optimization
**Objective:** Improve system scalability and computational efficiency.

- **Ensure Pipeline Supports Cloud Platform Operation:**
  - Deploy the pipeline to cloud platforms (such as **AWS** or **Google Cloud**) and use cloud services for large-scale data processing.
  - Configure cloud environments to ensure scalability of data processing and model training.
  - AWS (EC2, S3), Google Cloud, Docker, Kubernetes (containerization)

- **Complete Parallel Computing Optimization:**
  - Implement multi-threading or multi-process parallel computing to optimize data processing and model training efficiency, especially in federated learning environments.
  - Multiprocessing, Threading, Dask -->


### Phase 5: Automated Testing and Continuous Integration
**Objective:** Ensure code quality, automate testing, and integrate into the development process.

- **Set Up Automated Unit Testing and Continuous Integration Pipeline:**
  - Use **pytest** for unit testing to ensure the correctness of each functional module.
