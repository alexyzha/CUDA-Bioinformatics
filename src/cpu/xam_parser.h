#include "../sam_read.h"
#include "../util.h"

std::vector<sam_read*> read_sam(std::string file_path);

std::unordered_map<std::string, std::vector<std::string>> read_sam_header(std::string file_path);