#include "sam_container.h"
#include "util_structs.h"
#include "fx_util.h"

void sam_read_to_file(std::ofstream& out, sam_read& read);

void sam_to_file(sam_container& sam, std::string file_path);

std::string make_cigar(alignment& align);

std::vector<sam_read*> map_reads_to_ref(std::string ref, std::vector<fq_read*>& reads, size_t k);