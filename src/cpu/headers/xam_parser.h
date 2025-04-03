#include "sam_read.h"
#include "util.h"

/**
 *  @brief Reads all sam reads from a file.
 *  @param file_path `string`
 *  @return `vector<sam_read*>`
 *  @note Does not read headers.
 */
std::vector<sam_read*> read_sam(std::string file_path);

/**
 *  @brief Reads all sam headers from a file.
 *  @param file_path `string`
 *  @return `unordered_map<string, vector<string>>`
 *  @note Does not read sam reads.
 */
std::unordered_map<std::string, std::vector<std::string>> read_sam_headers(std::string file_path);