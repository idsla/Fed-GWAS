# Open issues

## 1.  globalseed data type is int getting 'numpy.int64' object has no attribute 'strip'. Below is the
elif self.current_stage == "init_chunks":
            try:
                global_seed_np = np.array([self.global_seed], dtype=np.int64)
                params = ndarrays_to_parameters([global_seed_np])
                self.current_stage = "iterative_king"
                return params, {} this is server code

partial_str = parameters[0].tobytes().decode("utf-8").strip()
        lines = partial_str.splitlines()
        logging.info(f"[Client {client.client_id}] Received {len(lines)} partial KING records from server.")
 
        for line in lines:
            # Expect: sampleA_ano sampleB_ano partial_phi n1_star
            parts = line.strip().split()
            if len(parts) < 4:
                continue this is code from client for iterative king                

## 2. two clients use the same file names, they rename or delete each other's files — race condition. Fixed it for init chunk need to fix it for exclude snps
If client folders are not isolated, there is conflict. 


def exclude_snps(plink_prefix, snp_list, new_prefix="filtered_data"):
    """
    Exclude the given SNPs from plink_prefix -> produce new_prefix .bed set.
    """
    if not snp_list:
        return plink_prefix
    exclude_file = "temp_exclude_snps.txt"

## 3. Error in KING records from server parameters: 'utf-8' codec can't decode byte 0xc9 in position 0: invalid continuation byte
