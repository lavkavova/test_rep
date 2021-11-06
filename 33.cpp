#include <iostream>
#include <cassert>

struct Node{
    int value;
    Node* next;
    Node(int v, Node* next = nullptr) : value(v), next(next){}
};

class LinkedList {
private:
    Node* head;
    int count;

public:
    LinkedList() : head(nullptr), count(0) {}
    ~LinkedList();

    void add(int v);
    void del();
    void remove(Node* after);
    void insert(Node* after, int v);
    Node* get_item(int index);
    int size() const { return count; }

    LinkedList(const LinkedList& source);
    void operator=(const LinkedList& source);

    void printList(){
        Node *current = head;
        while(current != nullptr){
            std::cout << current->value << " ";
            current = current->next;
        }
        std::cout << "\n";
    }
    void Ballon(int m, int k);
};

void LinkedList::add(int v) {
    if (head == nullptr){
        head = new Node(v);
    } else {
        head = new Node(v, head);
    }
    count++;
}

void LinkedList::del(){
    if (head == nullptr) { assert(false); }
    else if (count == 1){
        delete head;
        head = nullptr;
        count--;
    } else{
        Node* n = head->next;
        delete head;
        head = n;
        count--;
    }
}

void LinkedList::remove(Node *after) {
    Node* del_node = after->next;
    after->next = del_node->next;
    delete del_node;
    count--;
}

void LinkedList::insert(Node *after, int v) {
    Node* temp = new Node(v, after->next);
    after->next = temp;
    count++;

}

Node* LinkedList::get_item(int index) {
    int i = 0;
    Node* temp = head;
    while (temp != nullptr) {
        if (i == index)
            return (temp);
        i++;
        temp = temp->next;
    }
}

LinkedList::LinkedList(const LinkedList &source) {
    head = nullptr;
    count = 0;
    Node* n = source.head;
    LinkedList temp;
    while (n != nullptr){
        temp.add(n->value);
        n = n->next;
    }
    Node* t = temp.head;
    while (t != nullptr){
        add(t->value);
        t = t->next;
    }
}

void LinkedList::operator=(const LinkedList &source) {
    LinkedList n(source);
    std::swap(head, n.head);
    std::swap(count, n.count);
}

LinkedList::~LinkedList(){
    while(head != nullptr){
        Node* next_node = head->next;
        delete head;
        head = next_node;
    }
}

void LinkedList::Ballon(int m, int k){
    int j = 1;
    int n = size();
    for(int i = 0; i < k; ++i)
    {
        j = k;
        remove(get_item(j));
    }
}

int main() {
    LinkedList q;
    q.add(12);
    q.add(234);
    q.add(435);
    q.add(546);
    q.insert(q.get_item(0), 322);
//    q.remove(q.get_item(2));

    std::cout << q.size() << "\n";
    q.printList();
    LinkedList w = q;
    w.printList();
    std::cout << w.size() << "\n";
    // std::cout << w.get_item(2);
    std::cout << "\n";
    q.Ballon(3, 3);
    q.printList();

}