#ifndef LOGGER_H
#define LOGGER_H

#define LOG_MESSAGE_SIZE 256

void logger_info(const char *message);
void logger_error(const char *message);
void logger_debug(const char *message);

#endif /* LOG_H */