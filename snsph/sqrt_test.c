
float sqrtf_fast(float x);

main(int argc, char *argv)
{
	float y;
	float x = 100.0;

	y = sqrtf_fast(x);

	printf("answer is %f\n", y);

	exit(0);
}

