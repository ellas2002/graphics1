#include "FPToolkit.c"
#include <stdio.h>

void swap(int *a, int* b){
  int temp = *a;
  *a = *b;
  *b = temp;
}


void sort(int v[], int n) {
    // Selection sort
    for (int i = 0; i < n - 1; i++) {
        int min = i;
        for (int j = i + 1; j < n; j++) {
            if (v[j] < v[min]) {
                min = j;
            }
        }
        // Swap the found minimum element with the first element
        if (min != i) {
            swap(&v[min], &v[i]);
        }
    }
}

void printArray(int arr[], int size)
{
    int i;
    for (i = 0; i < size; i++){printf("%d ", arr[i]);}
    printf("\n");
}

int main(){
  int arr[] = {64, 25, 12, 22, 11};
  int n = sizeof(arr)/sizeof(arr[0]);
  sort(arr, n);
  printf("Sorted array: \n");
  printArray(arr, n);


  return 0;

}
