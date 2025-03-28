#include "sam_read.h"
#include "util.h"

/*
 *  Reads all sam reads from a file. Does not read headers.
 *  @param file_path `string`
 *  @return `vector<sam_read*>`
 */
std::vector<sam_read*> read_sam(std::string file_path);

/*
 *  Reads all sam headers from a file. Does not read sam reads.
 *  @param file_path `string`
 *  @return `unordered_map<string, vector<string>>`
 */
std::unordered_map<std::string, std::vector<std::string>> read_sam_headers(std::string file_path);