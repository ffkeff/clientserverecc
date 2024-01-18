#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <arpa/inet.h>

#define CLIENTS_COUNT 2
#define PORT 7777

void
listener(const int sfd, int cfd[2]){

	int max_sd, activity;
	fd_set readfds;

repeat:

	FD_ZERO(&readfds);
	FD_SET(sfd, &readfds);
	max_sd = sfd;

	for(int i = 0; i < CLIENTS_COUNT; i++){
		int sd = cfd[i];
		if(sd > 0)
			FD_SET(sd, &readfds);

		if(sd > max_sd)
			max_sd = sd;

	}

	activity = select(max_sd + 1, &readfds,
					  NULL, NULL, NULL);

	if(FD_ISSET(sfd, &readfds)){
		int new_socket;
		struct sockaddr_in clientaddr = {0};
		socklen_t clientaddrlen = sizeof(clientaddr);
		new_socket = accept(sfd,
							(struct sockaddr*) &clientaddr,
							&clientaddrlen);

		for(int i = 0; i < CLIENTS_COUNT; i++)
			if(cfd[i] == 0){
				cfd[i] = new_socket;
				break;
			}
	}

	for(int i = 0; i < CLIENTS_COUNT; i++){
		int sd = cfd[i];
		if(FD_ISSET(sd, &readfds)){
			char buffer[1024];
			if(read(sd, buffer, 1024) != 0){
				printf("[Send] Client [%d]: %s", i+1, buffer);
				int dest = (i + 1) % CLIENTS_COUNT;
				if(cfd[dest] != 0)
					write(cfd[dest],
						  buffer, strlen(buffer)+1);
			} else {
				close(sd);
				cfd[i] = 0;
			}
		}
	}

	goto repeat;
}

int
main(void){

	/*create and bind server socket*/
	int sfd = socket(AF_INET, SOCK_STREAM, 0);
	if(sfd < 0){
		fprintf(stderr, "[Error] socket()\n");
		return EXIT_FAILURE;
	}

	struct sockaddr_in servaddr = {
		.sin_family = AF_INET,
		.sin_port = htons(PORT)
	};

	socklen_t servaddrlen = sizeof(servaddr);
	if(bind(sfd, (struct sockaddr*) &servaddr, servaddrlen) < 0){
		fprintf(stderr, "[Error] bind()\n");
		return EXIT_FAILURE;
	}

	if(listen(sfd, CLIENTS_COUNT) < 0){
		fprintf(stderr, "[Error] listen()\n");
		return EXIT_FAILURE;
	}

	fprintf(stdout, "[OK] Server started\n");

	int cfd[CLIENTS_COUNT] = {0};
	listener(sfd, cfd);

	/*close sockets*/
	for(int i = 0; i != CLIENTS_COUNT, i++)
		if(cfd[i] == 0)
			close(cfd[i]);

    close(sfd);

    return 0;
}

