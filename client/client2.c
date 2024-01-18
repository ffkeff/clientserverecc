#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>
#include <assert.h>
#include <string.h>
#include "ecc.h"

#define HOST "185.185.68.109"
#define PORT "7777"

typedef struct{

	point_t *generator;
	mpz_t 	a;
	mpz_t 	b;
	mpz_t 	p;
	mpz_t 	cofactor;

} curve_t;

/*make socket by host and port*/
int
make_socket(char *host, char *port){

	struct addrinfo hints = {
		.ai_family = AF_INET,
		.ai_socktype = SOCK_STREAM,
		.ai_protocol = 0
	}, *serv_info, *p;

	int retval;
	if((retval = getaddrinfo(host, port, &hints, &serv_info)) < 0){
		fprintf(stderr, "getaddrinfo: %s\n", gai_strerror(retval));
		exit(EXIT_FAILURE);
	}

	int sock;
	for(p = serv_info; p != NULL; p = p->ai_next){
		if((sock = socket(p->ai_family, p->ai_socktype, p->ai_protocol)) < 0){
			close(sock);
			continue;
		}

		if(connect(sock, p->ai_addr, p->ai_addrlen) < 0){
			close(sock);
			continue;
		}
		break;
	}

	return sock;
}

point_t
*calc_point(curve_t curve, char *generator_x, char *generator_y, char *k){

	curve.generator = (point_t*)malloc(sizeof(point_t));

    const char *a = "6277101735386680763835789423207666416083908700390324961276";
    const char *b = "2455155546008943817740293915197451784769108058161191238065";
    const char *p = "6277101735386680763835789423207666416083908700390324961279";

	mpz_t secure_key;
    mpz_init(curve.generator->x); mpz_init(curve.generator->y);
    mpz_init(curve.a); mpz_init(curve.b); mpz_init(curve.p); mpz_init(curve.cofactor);
    mpz_init(secure_key);

	mpz_set_str(curve.generator->x, generator_x, 10);
    mpz_set_str(curve.generator->y, generator_y, 10);
    mpz_set_str(curve.a, a, 10);
    mpz_set_str(curve.b, b, 10), mpz_set_str(curve.p, p, 10);

	mpz_set_str(secure_key, k, 10);

	point_t *public_key = scalar_multiply(curve.generator, secure_key, curve.a, curve.b, curve.p);

	return public_key;
}

int
main(void){

	/*create connection with server*/
	int cfd = make_socket(HOST, PORT);
	fprintf(stdout, "[Success] connect to HOST\nStarting generation...\n");

	/*gen public key*/
	curve_t curve;
	char *key = "6";
	char *generator_x = "602046282375688656758213480587526111916698976636884684818";
    char *generator_y = "174050332293622031404857552280219410364023488927386650641";
	point_t *public_key = calc_point(curve, generator_x, generator_y, key);
	gmp_printf("public key: %Zd;%Zd", public_key->x, public_key->y);

    /*send public key to client*/
	char *pkey_x = mpz_get_str(NULL, 10, public_key->x);
	char *pkey_y = mpz_get_str(NULL, 10, public_key->y);
	char pkey[1024];
	strcpy(pkey, pkey_x);
	strcat(pkey, ";");
	strcat(pkey, pkey_y);
	strcat(pkey, ";");
	write(cfd, pkey, strlen(pkey) + 1);

    /*read public key from client*/
	char skey[1024];
	read(cfd, skey, strlen(skey) + 1);
	char *skey_x = strtok(skey, ";");
	char *skey_y = strtok(skey, ";");

	point_t *session_key = calc_point(curve, skey_x, skey_y, key);
	gmp_printf("session key: %Zd;%Zd", session_key->x, session_key->y);

	return 0;
}



