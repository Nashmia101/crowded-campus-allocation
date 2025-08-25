


# QUESTION 1

from collections import deque

class Flow_Vertex:
    """
    class description:
        
        The Flow vertex class represents a vertex in the flow network. As this flow network would be used 
        to later build a residual network. A vertex in this flow network represents either student, time slot
        class, source, sink, super source or super sink. Each vertex also maintains the list of outgoing edges 
        that helps in connecting it to other vertices which in return allows traversal and flow to pass through
        the network.
        
    Attributes:
        id (int)
        This is a unique identifier assigned to the vertex which helps identify the vertex.
        
        edges:
             These are the list of flow edges that originate from this particular vertex. These
             edges help connect this particular vertex to the other vertices in the flow network
    
    """
    def __init__(self, id):
        
        """
        Function description:
            Creates a flow vertex object with a specified id and then stores its outgoing edges. It does this by
            setting a vertex ID and intializing a empty list to store the outgoing edges from it.
        
        Input:
            id (int):
                This is a unique identifier assigned to the vertex which helps identify the vertex.
                
        Postcondition:
            it is only creating flow vertex object 
            
        Time complexity:
            O(1)
        Time complexity Analysis:
            The time complexity is O(1) because creating a empty list and intializing a variable is 
            a constant time operation
            
        Space complexity:
            O(1)
            
        Space complexity analysis
           O(1)
        Input Space : O(1) because id holds a single integer 
          
        Auxiliary space is constant because no addtional space is used as only a empty list is being created
        """
        self.id = id
        self.edges = []



class Flow_Edge:
    """
    The Flow edge class represents a directed edge in a flow network and each edge is created between the
    a particular vertex to the destination vertex (like from u->v). hence it is an outgoing edge that connects 
    the particular vertices. Each edge comes with its capacity and flow. Flow will be constantly changed when 
    the flow network undergoes ford fulkerson. Each edge also holds a reference to its connecting edge in the residual
    graph. This helps update the flow network when there are any changes in ternms of flow in residual
    graph as it allows the flow to be adjusted easily.
    
    """
    def __init__(self, v, capacity, flow):
        """
        function description:
            
            The Flow edge class represents a directed edge in a flow network and each edge is created between the
            a particular vertex to the destination vertex (like from u->v). hence it is an outgoing edge that connects 
            the particular vertices. Each edge comes with its capacity and flow. Flow will be constantly changed when 
            the flow network undergoes ford fulkerson. Each edge also holds a reference to its connecting edge in the residual
            graph. This helps update the flow network when there are any changes in ternms of flow in residual
            graph as it allows the flow to be adjusted easily. Hence allowing both the residual and flow graph to be linked.
            
        Input:
            v (Flow Vertex):
                This is the desitination vertex to which there is the outgoing edge. It is a Flow Vertex object
            capacity (int):
                This is the maximum flow that can be pushed through an edge. This value is fixed
            flow (int):
                This is the flow through that is going through the edge. This can be adjusted hence not fixed
                
        Postcondition:
            an Edge object is created with attributes v,capacity,flow and residual edge which helps link it to residual graph
            
        Time complexity:
            O(1)
            
        Time complexity analysis:
            Time complexity is O(1) because the constructor performs constant operation of assigning values
            
        :Space complexity: O(1)
       
        :Space complexity analysis:
            
          Input Space : O(1) because all inputs like v, capacity, flow  hold single values
          
          Auxiliary space is constant because no addtional space is used as the number of attributes are fixed

        """
        self.v = v
        self.capacity = capacity
        self.flow = flow
        self.residual_edge = None


class Flow_Graph:
    """
    class description:
        
        This class creates a Flow network with source, students, timeslots, classes and sink. This flow network helps
        us solve an allocation problem in which we have to allocate 2004 applied students to classes at different timeslots.
        This flow network connects source to students and students to timeslots and and timeslots to classes with lower bounds
        and maximum capacity and demands. Then it eliminates the lower bounds  and adjusts the demands and then connects the flow network 
        with super source and super sink. 
    """
    def __init__(self, n, m, timePreferences, proposedClasses):
        """
        Function description:
            This constructor creates a flow network by creating vertices and storing edges with lower bounds and max capacity
            The network is not fully constructed until we get rid of the demands and lower bouunds.
            What we first do is that we create a seperate adjacency list for every category like student, timeslot and classes and 
            store their vertex object and outing edges there. It records the edges connecting the vparticular vertices with lower bound
            and capacity. These edges are not yet finalized (connected) but are juststored so that we can process them later in a seperate 
            method that handles lower bounds and adjusts demands 
            
        Approach Description:
            The constructor first begins by creating vertex objects for each students, timeslots which is 
            fixed at 20 and classes and also for source. Then it stores these vertex objects for each category in 
            different adjacency lists along with their outgoing edges.
            Instead of creating Flow edge object immediately we first record the edges with the lower bound
            and capacity. These recorded edges are stored in the self.edges_with_lower_bounds
            edges from the source to the students are recorded with the lower bound of 1 and capacity of 1
            to ensure that every student must get assigned to atleast one class. edges from students to timeslots 
            are recorded to with lower bound of 0 and capacity of 1 so that students can have the option of choosing whatever 
            timeslot they want instead of restricting them. Timeslots to classes edges are recorded to have
            lower bound of minimum students that can be enroled to that class and capacity is the maximum students that can be allocated to 
            that class.
            There is also a demand array which is intialized to track demand for all vertices and the 
            source node is intialized the demand of -n so that it pushes flow to n number of students
            Each class node has a demand equal to the capacity and all other nodes are assigned 0.
            
        Input:
            n(int) is the number of students
            m(int) is the number of classes 
            timePreferences is the list of n lists where timePreferences[i] contains a permutation of the elements of set {0,1,...,19} to indicate the time slot preferences of student i
            proposedClasses is the list of m lists where each proposedClasses[i] represnets a list with [timeslot, minimum students to be allocated, maximum students to be allocated]
            
        Postcondition:
            flow network's vertices are created and stored 
            demand array is created
            edges are stored with data (lower bound and capacity)
            
        Time complexity:
             O(n + m)
             where n is the number of students and m is the number of classes
             
        Time complexity analysis:
            where n is the number of students and m is the number of classes
            
            O(n): when we create the student vertices and store them in self.student_nodes
            O(20): when we create timeslot vertices and store them in self.timeslot_nodes
            O(m): when we create class vertices and store them in class_nodes
            O(n x 20): creating edges from students to timeslots there are going to be 20 edges 
            because timeslots is bounded by 20
            O(m): when we create class we will create one edge to timeslot 
            O(n): when we add edge between source and student 
            O(n + m): when we initialize and assign demand values this is because we have allocate
            the demand arrazy size of (1+n+20+m)
               and then setting the source demand is a constant opertion of O(1)
               And then setting the class demands is O(m)
              Hence the total time complexity is O(n) + O(20) + O(m) + O(nx20)+ O(m) + O(n) + O(n + m)
              HEnce when we simplify this we get O(n + m)
              we will remove m from the time complexity because we are not allowed to make any assumptions 
              about whether m is greater then or smaller then n according to the assignment specification 
              hence will collapse  O(n + m) to O(n)
              
              
        Space complexity:
            
            O(n + m)
            where n is the number of students and m is the number of classes
            
        Input Space complexity:
            O(n) for the timePreferences list
            O(m) for the proposedClasses list
            
            Total Input space complexity is O(n + m)
            
        Auxiliary Space complexity:
            
            O(n): when we create a adjacency list of size n for students (self.student_nodes)
            O(1): when we create a adjacency list for timeslots that are bounded by
            20 hence this is a constant operation (self.timeslot_nodes)
            O(m): when we create a adjacency list of size m for classes (self.class_nodes)
            O(n x 20) when we store edges with lower bound and capacity from students to timeslots
            O(m): when we store edges with lower bound and capacity from timeslots to classes
            O(n): when we store edges with lower bound and capacity from source to student
            O(n + m): when we create a demand array of size (1+n+20+m)
            
            hence the total auxiliary space complexity is O(n + m)
       Total space complexity is O(n + m)
       we will remove m from the time complexity because we are not allowed to make any assumptions 
       about whether m is greater then or smaller then n according to the assignment specification 
       hence will collapse  O(n + m) to O(n)
            
        """
        self.n = n
        self.m = m
        self.timePreferences = timePreferences # intializes adjacency lists for for different categories
        self.proposedClasses = proposedClasses
        self.student_nodes = []
        self.timeslot_nodes = []
        self.class_nodes = []
        self.edges_with_lower_bounds = [] # edges with min and upper bound

        self.source_id = 0 # creates source id
        self.source_node = Flow_Vertex(self.source_id)
        self.super_source_id = -1
        self.super_sink_id = -1
        self.super_source_vertex = None
        self.super_sink_vertex = None

        for i in range(n): # creating student vertices and storing them in a adjacency list with their outgoing edges 
            v = Flow_Vertex(i + 1)
            self.student_nodes.append((v, []))

        for j in range(20):  # creating timeslot vertices and storing them in a adjacency list with their outgoing edges
            b = Flow_Vertex(n + j + 1)
            self.timeslot_nodes.append((b, []))

        for k in range(m):  # creating classes vertices and storing them in a adjacency list with their outgoing edges
            a = Flow_Vertex(n + 20 + k + 1)
            self.class_nodes.append((a, []))

        self.num_vertices = 1 + n + 20 + len(proposedClasses) # intail number of vertices 
        self.demand = [0] * self.num_vertices # Intialize demand array 
        self.demand[self.source_id] = -n #  intial demand

        # setting the demand of the class to its maximum capacity
        for class_id in range(len(proposedClasses)):
            class_vertex = proposedClasses[class_id]
            class_vertex_id = n + 1 + 20 + class_id
            self.demand[class_vertex_id] = proposedClasses[class_id][2]
            
        # Build student to timeslot edges
        for t in range(len(self.timePreferences)):
            time_prefs = self.timePreferences[t]
            student_vertex, adj_list = self.student_nodes[t]
            for p in range(len(time_prefs)):
                pref_timeslot = time_prefs[p]
                time_slot_vertex, c = self.timeslot_nodes[pref_timeslot]
                adj_list.append(time_slot_vertex)
                self.add_edge_with_lower_bounds(student_vertex.id, time_slot_vertex.id, 0, 1)

        # Build timeslot to class edges 
        for c_id in range(len(self.proposedClasses)):
            class_id = self.proposedClasses[c_id]
            future_class = class_id[0]   # timeslot index
            time_vertex, time_adj = self.timeslot_nodes[future_class]
            class_vertex, c = self.class_nodes[c_id]
            time_adj.append(class_vertex)
            lower = class_id[1]
            cap = class_id[2]
            self.add_edge_with_lower_bounds(time_vertex.id, class_vertex.id, lower, cap)

        # Build source to student edges
        for i in range(len(self.student_nodes)):
            student_vertice = self.student_nodes[i][0]
            self.add_edge_with_lower_bounds(self.source_id, student_vertice.id, 1, 1)

    def add_edge_with_lower_bounds(self, u, v, lower, cap):  # this stores lower bound and capacity for each edge 
        """
        Function description:
            Records an edge with between two vertices. this particular edge has the lower bound and the maximum
            capacity. Then we store it in a list for processing later on. This function does not create a
            list rather postpones creation of a list until the lower bound and demands resolved
        Inputs
         u(int): The id of the source vertex
         v(int): The id of the destination vertex
         lower(int): The lower bound 
         cap (int): The capacity
         
         Postcondition:
             The edge is stored in the list self.edges_with_lower_bounds
             
         Time complexity:
             O(1)
             appending a tuple to the list is a constant time operation which is O(1)
        Space complexity:
            Input space:
                u,v,lower,cap all are single value integers hence the input space will be O(1)
            Auxiliary Space:
                no additional storage used so O(1) auxiliary space
        
        """
        self.edges_with_lower_bounds.append((u, v, lower, cap))

    def eliminate_lower_bound_and_adjust_demand(self): 
        """
        Function description:
            converts the flow network by eliminating lower bounds on edges and adjusting the vetices demands
            This is important in order to for converting the circulation with the lower bounds into the standard flow network
            for each edge that goes from a source vertex (u) yo destination vertex (v)
            we subtract lower bound from the demand of u and then add the lower bound to the demand of v
            and then we replace the original edge with the capacity of capacity - lower bound. 
            After doing this we would connect a super sourc node is connected to all the vertices 
            with negative net demand
            and we would connect a super sink node to the classes with the positive net demand 
            
            Postcondition:
                All lower bounds are removed from the graph 
                The graph contains only regular edges with capacity
                the demand array contains the the inflow/outflow for each vertex
                the super source and super sink are added in the graph now 
                
        Input:
            None
            
        Time complexity:
            
            O(n + m)
            where n is the number of students and m is the number of classes
            
        Time complexity Analysis:
            where n is the number of students and m is the number of classes
            
            - O(n + m): Intializing the fix_bound array with the size (1+n+m+20) hence it will be 
                      O(n+m)
            - O(n + m):
                    then we will iterate over the edges self.edges_with_lower_bounds in order to resolve the 
                    the lower bounds 
                    which will lead to O(E) complexity where E is the number of edges 
                    E includes upto (n x 20) edges between student and timeslots as timeslots are
                    bounded by 20 and it will also include upto n edges from source to students
                    It also includes edges from timeslot to classes which can be upto m number of edges
                    Hence we can say that O(E) is the O(n + m) ( as 20 is bounded and constant)
                    
           - O(n + m): we'll also adjust the demand array. This will be O(n + m) because we'll
           have to update the demand of each vertex 
           
           - O(n + m): when we'll create a self.edges_with_bounds list which will lead to O(n + m)
                       because we will add each edge which will take O(n + m)
           - O(m):     to create a demand_class boolean array in order to know where to connect the super sink 
           - O(n + m): when we'll connect the super source and super sink  according to the demand.
                       It will be O(n + m) because we'll have to check all the vertices to see there demand
           - O(n + m) : when we'll call the function add_edge because its time complexity is O(n + m)
           
           Hence the total time complexity is O(n + m)
           we will remove m from the time complexity because we are not allowed to make any assumptions 
           about whether m is greater then or smaller then n according to the assignment specification 
           hence will collapse  O(n + m) to O(n)
           
        Space complexity
           
            O(n + m)
            where n is the number of students and m is the number of classes
            
        Space complexity analysis:
            
            O(n + m): fix_bound arrayhas been created with the size num_vertices which is (n + m) because 
                      20 is bounded and there is only one source vertices which gives us (n + m)
            O(n + m): edges_with_bounds list is created to store the regular edges after removing lower 
                      bounds. They are the same size as edges_with_lower_bounds O(n + m)
            O(n + m): boolean array created with the size of n + m + 20 + 1 so thus it will be O(n + m) because rest are constant
            O(n + m): add_edge method" space complexity
            
        So total space complexity is  O(n + m)
        we will remove m from the time complexity because we are not allowed to make any assumptions 
        about whether m is greater then or smaller then n according to the assignment specification 
        hence will collapse  O(n + m) to O(n)
        
        """
        fix_bound = [0] * self.num_vertices  # intializing a list to store the net demand at each node after elimination of the lower bound

        for u, v, lower, cap in self.edges_with_lower_bounds: # adjust demands according to the lower bounds
            fix_bound[v] = fix_bound[v] + lower
            fix_bound[u] = fix_bound[u] - lower

        for v in range(self.num_vertices): # subtract the fixed demand from each node
            self.demand[v] = self.demand[v] - fix_bound[v]

        # create super and sink vertices 
        self.super_source_id = self.num_vertices
        self.super_sink_id = self.num_vertices + 1
        self.super_source_vertex = Flow_Vertex(self.super_source_id)
        self.super_sink_vertex = Flow_Vertex(self.super_sink_id)

        self.edges_with_bounds = []  
        for u, v, lower, cap in self.edges_with_lower_bounds: # convert the lower bound edges into the regular edges after max - min
            remain = cap - lower
            self.edges_with_bounds.append((u, v, remain))
        
        demand_class = [0] * self.num_vertices
        for class_node, node in self.class_nodes:
            demand_class[class_node.id] = 1
        
        # Step 2: Only check demand > 0 and if node is a class
        for v in range(self.num_vertices):
            d = self.demand[v]
            if d < 0:
                self.add_edge(self.super_source_id, v, -d, 0)
            elif d > 0 and demand_class[v]:
                self.add_edge(v, self.super_sink_id, d, 0)
                        

        for u, v, cap in self.edges_with_bounds:
            self.add_edge(u, v, cap, 0) # add the edges after eliminating lower bound and adjusting demands 

    def add_edge(self, u, v, capacity, flow): # connect edges between two vertices 
        """
        Function description 
         This method connects two vertices in the flow network by connecting an edge between the 
         two vertices with the specified capacity and flow. It first gets Flow vertex object 
         for the 2 vertices  and then it appends the Flow_edge object into to the adjacency list (edges)
         of the first vertex.
         
         Input
         
         u(int): id of the source vertex of from which the edge comes
         v(int): id of the destination vertex to which the edge points
         capacity(int): The maximum possible flow that can be pushed through this edge
         flow(int): The current flow assigned to the edge
         
         Postcondition:
            creates edge between vertices 
            
        Time complexity:
            O(n + m)
            where n is the number of students and m is the number of classes
        Time complexity Analysis:
            O(m + n): when the get_vertex function is called because its time complexity is O(n+m)
            O(1): a new Flow object is created which is a O(1) operation
            o(1): appending the object to the adjacency list of the particular vertex is also a 
                  constant time operation
          we will remove m from the time complexity because we are not allowed to make any assumptions 
          about whether m is greater then or smaller then n according to the assignment specification 
          hence will collapse  O(n + m) to O(n)
        
        Space complexity:
              O(1)
        Space complexity Analysis:
            
            Input space: u, v, capacity, flow are single interger values hence they take up no additional space so it is a constant O(1)
            Auxiliary Space: No additional storage is being used so O(1) space
        """
        first_vertex = self.get_vertex(u)
        second_vertex = self.get_vertex(v)
        if first_vertex is not None and second_vertex is not None:
            first_vertex.edges.append(Flow_Edge(second_vertex, capacity, flow))

    def get_vertex(self, id):
        """
        Function description: 
            This function returns the Vertex object at a particular index in the all_vertices list
        
        :Input:
            
            id: int ( it is the index of the vertex in the all_vertices list)
            
             
        :return:
           
            returns the appropriate Vertex object
         
        Time complexity:
            O(n + m)
            where n is the number of students and m is the number of classes
            
        Time complexity analysis:
            O(1): all_vertices[id] is a constant time operation
            
            O(n + m): then building the all_vertices array takes up O(n + m) because you go append
                      all student vertices n, 20 timeslots, source vertex and class vertices m and super 
                      source and super sink if they exist. This builds upto O(n + m) because we have to iterate over n,m,20 so it will
                      be O(n + m) because rest is constant
          Total time complexity is O(n + m)
          we will remove m from the time complexity because we are not allowed to make any assumptions 
          about whether m is greater then or smaller then n according to the assignment specification 
          hence will collapse  O(n + m) to O(n)
          
          Space complexity:
              Input space:
                  Input Space : O(1) because id holds a single integer
              Auxiliary Space:
                  O(n + m)
              where n is the number of students and m is the number of classes
              then building the all_vertices array takes up O(n + m) because you go append
              all student vertices n, 20 timeslots, source vertex and class vertices m and super 
              source and super sink if they exist. This builds upto O(n + m) because rest is constant 
         Total space complexity is O(n + m)
         we will remove m from the time complexity because we are not allowed to make any assumptions 
         about whether m is greater then or smaller then n according to the assignment specification 
         hence will collapse  O(n + m) to O(n)
        """
        all_vertices = [self.source_node] # storing every vertex after here the elimination and adjustment process of demands 
        for su, st in self.student_nodes:
            all_vertices.append(su)
        for ti,tim in self.timeslot_nodes:
            all_vertices.append(ti)
        for cl,cla in self.class_nodes:
            all_vertices.append(cl)
        if self.super_source_vertex is not None and self.super_sink_vertex is not None:
            all_vertices.append(self.super_source_vertex)
            all_vertices.append(self.super_sink_vertex)
            
        return all_vertices[id]
            


class Res_Vertex:
    """
    Representa a vertex in the residual network. This residual network is used in Ford Fulkerson
    """
    def __init__(self, id):
        """
        Function description:
            Constructs a residual graph vertex and intializes its attributes 
            Attributes that are intialized 
            
        Input:
            id: int ( it is the index of the vertex in the all_vertices list)
            
        Postcondition:
            a Res_Vertex vertex object is created  with the given attributes
        
        Time complexity:
            O(1)
        Time complexity analysis:
            Time complexity is O(1) because the constructor performs constant operation of assigning values
        Space complexity:
            O(1)
        Space complexity analysis:
            Input Space : O(1) because id holds a single integer
            
            Auxiliary Space: Auxiliary space is constant because no addtional space is used as the number of attributes are fixed
                    
        
        """
        self.id = id
        self.edges = []      # List of Redge
        self.prev = None     # Used for backtracking
        self.visited = False
        self.min_flow = float('inf')


class Redge:
    """
    class description
        Represents a residual edge in a residual graph used for flow augmentation in ford fulkerson algorithim
    
    """
    def __init__(self, u, v, for_weight, back_weight):
        """
        Function description:
            This function intializes a residal edge object between two vertices in the residual graph
            it also stores backward and forward residual capacities. It also maintains a referernes to the original 
            Flow_Edge from the flow graph allowing updates to be made to the actual flow values when the residual 
            graph is updated
        
        Input:
            u (Res_vertex): the source vertex
            v (Res_vertex): the destination vertex
            for_weight (int): forward residual capacity from u to v
            back_weight (int): Backward residual capacity from v to u
            
        Time complexity:
            O(1)
            
        Time complexity Analysis:
            Time complexity is O(1) because the constructor performs constant operation of assigning
            values
        
        Space complexity:
            O(1)
        
        Space complexity Analysis:
          Input Space : O(1) because all inputs like u, v, for_weight, back_weight hold single values
          
          Auxiliary space is constant because there are single value integers and list of length 2

        """
        self.v = v
        self.u = u
        self.w = [for_weight, back_weight]
        self.flow_edge = None # points to the flow edge in the flow network 

class Residual:
    """
    Function Description:
        This constructor builds the residual graph from the given flow network
        Which is then used in the ford fulkerson to augmnet the path which guarantees 
        that every constraint is fulfilled and every student gets allocated to a class. 
        For each edge in the flow network it creates two coressponding edges a 
        forward residual edge with capacity ( capacity - minimum flow) and a 
        backward residual edge with capacity equal to the current flow
        These residual edges allow the Ford-Fulkerson algorithm to find augmenting paths and push flow
        through them They can even reverse the flow when needed.
        
    Approach Description:
        First it creates th residul graph by creating residual vertex(Res_Vertex) for each 
        vertex in the flow network and also creates residual vertex for super source and super sink.
        Then it iterates over all the vertices in order to iterate over th edges of those vertices
        then for each edge it calculates the forward residual capacity 
        and backward residual capacity and then adds both these to the residual network in order 
        to enable tp push the flow forward and backward
        
    Input:
        flow_graph:
            The original graph that contains every vertex  and edges with flow and capacity
            
    Postcondition
        Creates a Residual graph of the Flow network and adds backward and forward edge for 
        each edge in the flow network
        Enables forward and backward traversal of flow for augmenting path search
        
    Time complexity:
        O(n + m)
        where n is the number of students and m is the number of classes
        
    Time complexity analysis:
        Since we iterate over all the vertices in the original flow network in order to create 
        Res_Vertex object for each of those vertices. which will lead to O(n + m) as n vertex objects 
        will be created for the students and m vertex objects will be created for classes
        and 20 vertex objects will be created for the timeslots and then will create vertex objects for 
        source and super source and super sink so hence O(n + m + 20 + 1 +1+1) it will be O(n + m) because timelsots are bounded by 20
        Then the same way we'll iterate over edges for each vertex which will lead to O(E) complexity
        where E is the total number of Edges ( as we have edges from super source to source and then
        from source to student we'll have (n) edges and then we'll have 20 edges from student to timeslots because timelots are bounded
        by 20 and then from timeslots to classes we'll have (m) edges and then finally from classes to sink we'll have
        m edges was but since Timeslots is O(n + m + 20n + m + 1) hence we'll have O(m + n)
        Each edge results in 2 calls to add_edge() one for forward edge and one for one backward edge.
        Although add_edge() iterates through vertex's outgoing edge list (O(E_u + E_v)) the total
        work done across all vertices is still O(E) which is O(n + m) as discussed previously because the sum of all vertex degrees
        is bounded by the total number of edges hence the overall time remains O(n + m)
        
        Hence total time complexity is O(n + m)
        we will remove m from the time complexity because we are not allowed to make any assumptions 
        about whether m is greater then or smaller then n according to the assignment specification 
        hence will collapse  O(n + m) to O(n)
        
        
    Space complexity 
        O(n + m)
        where n is the number of students and m is the number of classes
    Space complexity analysis:
        
        Input Space:
             O(n + m) for the flow_graph
        Auxiliary Space:
            storing each and every vertex in the self.vertices will lead to a O(n + m) 
            auxiliary space
            
    so the total space complexity is O(n + m)
    we will remove m from the time complexity because we are not allowed to make any assumptions 
    about whether m is greater then or smaller then n according to the assignment specification 
    hence will collapse  O(n + m) to O(n)
    """
    def __init__(self, flow_graph):
        self.flow_graph = flow_graph
        self.vertices = [] # Initialize empty list to store Res_Vertex
        for i in range(flow_graph.super_sink_id + 1):  # adding Res_Vertex objects to the list
            self.vertices.append(Res_Vertex(i))
        for u in range(flow_graph.super_sink_id + 1):
            vertex = flow_graph.get_vertex(u)
            if vertex is not None:
                for edge in vertex.edges: # For each outgoing edge from this vertex
                    residual_capacity = edge.capacity - edge.flow # Forward edge residual capacity
                    flow_amount = edge.flow # Backward edge residual capacity

                    # Add both forward and backward residual edges
                    self.add_edge(u, edge.v.id, residual_capacity)
                    self.add_edge(edge.v.id, u, flow_amount)

    def get_vertex(self, id):
        """
        Function description: 
            This function returns the Vertex object at a particular index in the self.vertices list
        
        :Input:
            
            id: int ( it is the index of the vertex in the self.vertices list)
            
             
        :return:
           
            returns the appropriate Vertex object
            
        :Time complexity:
            O(1)
            
        :Time complexity analysis:
            Time complexity is O(1) because it access through index is a constant operation 
            
        :Space complexity: O(1)
       
        :Space complexity analysis:
            
         Input Space : O(1) because id holds a single integer
          
         Auxiliary space is constant because no addtional space is used
     
        """
        return self.vertices[id]

    def add_edge(self, u_id, v_id, forward_capacity):
        """
        Function description:
            This function adds or updates the the residual edges between two vertices u_id, v_id
            in the residual graph. This includes both the forward edge
           (with residual capacity) and the backward edge (with reverse capacity) Also
           links the forward residual edge to the original flow edge if it exists this way with changes made
           to the residula graph can be automatically made to the flow in the flow network.
           
         Input:
             u_id:
                 the id (unique identifier of the source vertex)
            v_id
                the id (unique identifier of the destination vertex)
            forward capacity(int):
                the foward residual capacity from u to v
                
        Postcondition:
            Updates the edges list of both u and v in the residual graph.
            also links the residual edge to its corresponding flow edge in the original graph
        
        Time complexity:
            
            O(E_u + E_v)
            
            where E_u is the number of outgoing edges from u vertex 
            where E_v is the number of outgoing edges from v vertex 
            
        Time complexity analysis:
            First we will search for an existing forward edge from u to v hence we will
            be iterating over outgoing edges of u hence the complexity here will be O(E_u)
            Then we will search for an existing forward edge from v to u hence we will
            be iterating over outgoing edges of v hence the complexity here will be O(E_v)
            then in order to check if the forward edge in the residual graph matches the one in 
            flow graph then we'll to iterate over outgoing edges of u hence the complexity here will be O(E_u)
            hence the total time complexity is O(E_u + E_v + E_v) which simplifies to O(E_u + E_v)
            
        Space complexity:
            O(1)
            
        Space complexity analysis:
            
            No additional spsace is used except for creating an Redge object for the edges which is a constant time operation
            
            
        """
        
        u = self.get_vertex(u_id)
        v = self.get_vertex(v_id)

       # Search if a forward edge from u to v already exists
        forward_found = None
        for e in u.edges:
            if e.v.id == v_id:
                forward_found = e
                break

        if forward_found is not None:
            forward_found.w[0] = forward_capacity # If found update its forward capacity
        else:
            forward_found = Redge(u, v, forward_capacity, 0)
            u.edges.append(forward_found)   # If not found create a new forward edge

            flow_vertex = self.flow_graph.get_vertex(u_id) # Link to the original flow edge if it exists
            if flow_vertex is not None:
                for f in flow_vertex.edges:
                    if f.v.id == v_id:
                        forward_found.flow_edge = f   # referencing so that updates in the residual graph can be reflected in the flow network's flow
                        f.residual_edge = forward_found
                        break

        # Search if a backward edge from v to u already exists
        backward_found = None
        for e in v.edges:
            if e.v.id == u_id:
                backward_found = e
                break

        if backward_found is not None:
            backward_found.w[1] = forward_capacity # If found update its backward capacity
        else:
            v.edges.append(Redge(v, u, 0, forward_capacity)) # If not found create a new backward edge


def bfs(created_residual_network, source_id, sink_id):
    """
    Function Description:
        This function performs a Breadth-First Search (BFS) on the residual graph to find
        an augmenting path from the source to the sink. This is used in each iteration
        of the Ford-Fulkerson algorithm to determine if more flow can be pushed from the source
        to the sink.
        It updates the prev pointer of each visited vertex to reconstruct the path which helps us in back tracking
        and it also tracks the minimum flow capacity available along the path
        using the min_flow attribute of each vertex.

    Approach Description:
        - The function starts by resetting all vertices in the residual graph to an unvisited state,
          with no previous edge and infinite min_flow this helps start the bfs from start. 
          This ensures no leftover data from previous searches affects the current search
        - Then we start the BFS from the source by marking super_source vertex and placed into the BFS queue. 
          This queue is used to explore the graph level by level ensuring shortest augmenting
             paths which means that paths are first found according to the number of edges they use like for example a path with less number 
             number of edges is found first. This queue is ued to store all the vertices that have been discovered throughout the BFS traversal
        - While the queue is not empty the algorithm dequeues the current vertex u.
           - For every outgoing residual edge from u it checks the forward residual capacity
             to see if more flow can be pushed.
           - If the destination vertex v is unvisited and the residual capacity is positive:
             - Mark v as visited.
             - Record it as the previous edge in order to backtrack later.
             - Update min_flow to track the minimum residual capacity along the path from super source to that vertex
             - If the sink is reached during BFS we immediately return True indicating a path exists
             - If the BFS completes without reaching the sink, return False (no more paths).

    Inputs:
        created_residual_network (Residual):
            The residual graph built from the original flow graph using forward and backward residual capacities.

        source_id (int):
            The ID of the super_source vertex in the residual graph.

        sink_id (int):
            The ID of the sink vertex in the residual graph.

    Returns:
        Returns True if an augmenting path from super_source to sink exists otherwise returns False

    Time Complexity:
        O(n + m) where n is the number of students and m is the number of classes
        

    Time Complexity Analysis:
        - Resetting all vertex attributes (`visited`, `prev`, `min_flow`) takes O(n + m) time,
          because the number of vertices is proportional to students + classes + timeslots + super nodes.
        - The BFS loop visits each vertex once and checks each outgoing edge once.
        - The while loop for the queue take O(n+m) time because the number of vertices are
            (1 + n + 20 + m + 1 +1) which takes O(n+m) time because timeslots and source and super sink and super source are bounded by constant.
        - Popleft form the queue take O(1) time
        - The for loop for iterating through all edges takes
          - super_source to source O(1)
          - source to students: O(n)
          - students to timeslots: ≤ 20 per students so 20* O(n)
          - timeslots to classes: O(m)
          - classes to sink: O(m) 
           adding all these together simplifies to O(n+m) time
        - Since the total number of residual edges in the graph is also bounded by O(n + m),
          the BFS completes in O(n + m) time.
          we will remove m from the time complexity because we are not allowed to make any assumptions 
          about whether m is greater then or smaller then n according to the assignment specification 
          hence will collapse  O(n + m) to O(n)

    Space Complexity:
        O(n + m) where n is the number of students and m is the number of classes

    Space Complexity Analysis:
        - Input space is O(n + m) for the residual graph.
        - Auxiliary space is O(n + m) due to the `discovered` queue as it can store upto all the vertices in the graph.
    The total space complexity is O(n + m)
    we will remove m from the time complexity because we are not allowed to make any assumptions 
    about whether m is greater then or smaller then n according to the assignment specification 
    hence will collapse  O(n + m) to O(n)
    
   """
    # Resetting each vertex in the residual graph 
    for ve in created_residual_network.vertices:
        ve.prev = None
        ve.visited = False
        ve.min_flow = float('inf')

    discovered = deque() # Mark source vertex as visited and add it to the queue 
    source = created_residual_network.get_vertex(source_id)
    source.visited = True
    discovered.append(source)

    while len(discovered) > 0: #staring BFS
        u = discovered.popleft() # Dequeue a vertex 
        for edge in u.edges:
            v = edge.v
            if not v.visited and edge.w[0] > 0: # only check vertices if they are unvisted and if the edge is greater then 0
                v.visited = True
                v.prev = edge
                v.min_flow = min(u.min_flow, edge.w[0]) # Track minimum flow along this path
                if v.id == sink_id:  # If sink is reached we found a valid augmenting path
                    return True # if no path found return true
                discovered.append(v) # add to the queue 
    return False # if no path found return false 


def ford_fulkerson(ford_flow_graph):
    """
    Function Description:
        This function implements the Ford-Fulkerson algorithm to compute the maximum flow 
        in a flow network where students are assigned to classes under capacity constraints.
        It repeatedly searches for augmenting paths in the residual graph using BFS 
        and updates flows along these paths until no more augmenting paths exist.
        The function then checks whether the total flow satisfies all class demand requirements (the total flow is greater than or equal to the total demand). 
        If it does the updated flow graph is returned otherwise none is returned 

    Approach Description:
        - Firstly a Residual graph is constructed from the given flow graph. This includes forward 
          edges (for remaining capacity) and backward edges (for current flow that can be undone).
        - The algorithm initializes the maximum_flow to 0 and repeatedly calls BFS to find 
          augmenting paths from the super source to the super sink.
        - If a path is found the minimum residual capacity along the path is 
          determined by min_flow.
        - The flow is then augmented by this minimum residual capacity:
            - The forward edge's residual capacity is decreased.
            - The backward edge's capacity is increased.
            - The corresponding flow in the original flow graph is updated.
        - The minimum residual capacity value is added to `maximum_flow`.
        - This continues until no more augmenting paths exist.
        - Finally, the sum of all positive demands in the flow graph is compared with the 
          computed maximum_flow . If they match, the flow is feasible and the updated 
          flow graph is returned. If not, the function returns None.

    Input:
        ford_flow_graph (Flow_Graph):
            The flow graph consisting of students, timeslots, classes, and super source/sink.
            Edges contain capacity, flow.

    Output:
        Returns the modified flow graph if a valid allocation is found thus the graph is considered feasible however if no allocation is found it returns     None.

    Time Complexity:
        O(n × (n + m)) in the worst case 
        where n = number of students
        where m = number of classes

    Time Complexity Analysis:
        - Each iteration of the while loop runs a BFS to find an augmenting path. 
          The BFS runs in O(N + M) time because:
            - The number of vertices is bounded by n (students) + m (classes) + 20 (timeslots) + 2 (super nodes), 
              which simplifies to O(n + m) since 20 and 3 are constants.
            - The number of edges is also O(n + m) as each student connects to upto 20 timeslots and 
              each class is connected to one timeslot and the sink.
        - Each BFS iteration performs updates on the flow path in O(n) time (bounded by the length of the path).
        - The total number of iterations is at most O(n) and  the max flow is likely equal
                         to the number of students and the total required flow is at most the number of students.
        - Therefore, total time complexity is O(n × (n + m)).
        we will remove m from the time complexity because we are not allowed to make any assumptions 
        about whether m is greater then or smaller then n according to the assignment specification 
        hence will collapse O(n^2)

    Space Complexity:
        O(n + m)

    Space Complexity Analysis:
        - Input space: the flow graph occupies O(n + m) space.
        - Auxiliary space: residual graph stores one forward and one backward edge for each original edge along with (n + m) residual vertices,
          resulting in O(n + m) extra space.
        - BFS uses a queue and attributes like (visited, prev, min_flow) which take O(n + m) space.
     we will remove m from the time complexity because we are not allowed to make any assumptions 
     about whether m is greater then or smaller then n according to the assignment specification 
     hence will collapse O(n+ m) to O(n)
    """
    # creating residual graph from the original flow graph
    flow_residual_graph = Residual(ford_flow_graph)

    maximum_flow = 0 # Initializing the total flow in the network 
    source_id = ford_flow_graph.super_source_id
    sink_id = ford_flow_graph.super_sink_id
 
    while bfs(flow_residual_graph, source_id, sink_id): # While augmenting paths exist from super source to sink in residual graph
        flow_in_path = flow_residual_graph.get_vertex(sink_id).min_flow
        v = flow_residual_graph.get_vertex(sink_id)
        while v.prev is not None: # start from the sink and backtrack to super source
            edge = v.prev
            edge.w[0] = edge.w[0] - flow_in_path # Reduce the forward edge capacity by the flow we’re adding
            for re_edge in edge.v.edges:
                if re_edge.v.id == edge.u.id:
                    re_edge.w[0] = re_edge.w[0] + flow_in_path # Increasing the backward edge capacity
                    break
            if edge.flow_edge is not None: # If this edge is part of the original flow graph update its flow
                edge.flow_edge.flow = edge.flow_edge.flow + flow_in_path
            v = flow_residual_graph.get_vertex(edge.u.id)
        maximum_flow = maximum_flow + flow_in_path #Adding flow for this path to the maximum flow

    # After aigmenting paths finish we check if flow satisfies demand
    sum_demand = 0
    for d in ford_flow_graph.demand:
        if d > 0:
            sum_demand = sum_demand + d

    if maximum_flow < sum_demand: # If we could not meet all demands return None
        return None
    return ford_flow_graph # otherwise return the flow graph with the updated flow

def pre_assign_students(n, m, timePreferences, proposedClasses):
    """
    Function Description:
        This function performs a preprocessing step for the crowded campus allocation problem.
        It attempts to assign students to classes based on their top-5 time slot preferences in a way that helps 
        each class meet its minimum capacity requirement before running the main ford Fulkerson algorithm.
        The function returns:
          - a list of pre-assignments (student, timeslot, class),
          - a flag array indicating whether each student has been preprocessed or not
          - and the total number of students who received satisfactory assignments.

    Approach Description:
        - Initializes a mapping of each timeslot to the list of classes offered during that timeslot.
        - Iterates through each student to determine whether they can be pre-assigned to any class based on:
            - Whether a class exists in their top-5 preferred time slots.
            - Whether that class still requires more students to meet its minimum capacity.
        - If a match is found, the student is added to the list of pre-processed students and tracking variables 
          are updated accordingly to avoid reassigning the same student or overfilling a class.
        - This preprocessing ensures that At least minimumSatisfaction students get allocated to classes with class times that are
          within their top 5 preferred time slots

    Inputs:
        n (int): Number of students.
        m (int): Number of proposed classes.
        timePreferences is the list of n lists where timePreferences[i] contains a permutation of the elements of set {0,1,...,19} to indicate the time slot      preferences of student i
       proposedClasses is the list of m lists where each proposedClasses[i] represnets a list with [timeslot, minimum students to be allocated, maximum students to be allocated]

    Returns:
        Tuple:
            - pre_processed_students (list of tuples): (student_id, timeslot, class_id) for each pre-assignment.
            - preprocessed (list of int): Binary flags indicating if a student was preprocessed (1) or not (0).
            - num_satisfied (int): Count of students successfully assigned in preprocessing.

    Time Complexity:
        O(n × m)
    where n is the number of students and m is the number of classes

    Time Complexity Analysis:
        - Building timeslot_connect_class list takes O(m) as each class is processed once and assigned to a particular timeslot.
        - For each student which are denoted by n the function checks their top 5 preferred time slot which is constant hence it would be O(n)
        - For each preferred time slot, it accesses the list of classes scheduled in that slot.
          Hence total time O(m + (n × m)) which simplifies to (n × m)
          we will remove m from the time complexity because we are not allowed to make any assumptions 
          about whether m is greater then or smaller then n according to the assignment specification 
          hence will collapse O(m x n) to O(n^2)
    Space Complexity:
        O(n + m)
    where n is the number of students and m is the number of classes
    Space complexity Analysis:
    Input Space complexity:

    O(n) for the timePreferences list
            O(m) for the proposedClasses list
            
    Total Input space complexity is O(n + m)

    Auxiliary Space Complexity:
    - pre_processed takes upto O(n) space 
    - class_preprocessed_count takes up to O(m) space
    - timeslot_connect_class` O(m) since total number of classes across timeslots is m
    - pre_processed_students`: Up to O(n) assignments stored

     So the total space complexity is O(m) + O(m) + O(n) + O(n) which simplifies to O(m + n)
     Hence total space complexity is  O(m + n)
     we will remove m from the time complexity because we are not allowed to make any assumptions 
     about whether m is greater then or smaller then n according to the assignment specification 
     hence will collapse O(m + n) to O(n)
     
   """
    pre_processed_students = []
    preprocessed = [0] * n
    num_satisfied = 0
    class_preprocessed_count = [0] * m

    timeslot_connect_class = []
    for t in range(20):
        timeslot_connect_class.append([])
    for cl in range(m):
        timeslot = proposedClasses[cl][0]
        timeslot_connect_class[timeslot].append(cl)

    for stu in range(n): # Iterate over each student to try assigning them based on their top 5 preferences
        if preprocessed[stu]:
            continue # skip if this student is already preprocessed
        time_pref = timePreferences[stu]
        i = 0
        while i < 5 and i < len(time_pref): # index to iterate through top 5 preferences
            time_p = time_pref[i]
            for ct in timeslot_connect_class[time_p]:
                min_sat = proposedClasses[ct][1] # get the minimum required students for class ct
                if class_preprocessed_count[ct] < min_sat:
                    # If this class still needs students then assign this student
                    pre_processed_students.append((stu, time_p, ct))
                    preprocessed[stu] = 1 # mark student as assigned
                    class_preprocessed_count[ct] = class_preprocessed_count[ct] + 1
                    num_satisfied = num_satisfied + 1 # increment number of satisfied (assigned) students
                    break
            if preprocessed[stu]: # break the loop if student was already assigned
                break
            i =  i + 1

    return pre_processed_students, preprocessed, num_satisfied


def collect_remaining_students(n, preprocessed, timePreferences):
    """
    Function description:
        This function identifies students who were not pre-assigned during the preprocessing phase and collects
        both their indices and time preference lists for use in the second-phase flow allocation.

    Approach description:
        The function iterates through all n students. For each student, it checks whether they were assigned in
        the pre-assignment phase using the preprocessed array which stores 0 or 1 to indicate whether a particular student has already been assigned or not      if a student has been assigned then its 1 otherwise 0. If a student has not been preprocessed 
        their index is added to the remaining_students list and their corresponding
        time preferences are added to decreased_Preferences_time. This separation allows the second phase of the
        algorithm to focus only on unassigned students while maintaining the original time preference data for 
        constructing the flow graph.

    :Input:
        n (int): Total number of students .
        preprocessed (list): A list of size n indicating whether each student has already been assigned (1) or not (0).
        timePreferences is the list of n lists where timePreferences[i] contains a permutation of the elements of set {0,1,...,19} to indicate the time slot preferences of student i

    :Return:
        Returns a tuple of two lists:
            - remaining_students (list) this list contains the indices of students who are still unassigned.
            - decreased_Preferences_time (list of list): Their corresponding time slot preference lists.

    :Time complexity:
        O(n) where n is the number of students 

    :Time complexity analysis:
        - A single loop runs iterating through each student hence it is O(n)
        - For each student we check if they are assigned or not.
        - If the student is unassigned its index is stored in remaining_students and its timeslot preferences are stored in remaining_students These are
          are constant operations hence O(1)
        Hence the total time complexity is O(n)

    :Space complexity:
        O(n) where n is the number of students 
    
    Space complexity analysis:
      
    Input Space complexity:
         n: this is just the the total number of students hence it takes constant space 
         preprocessed: this is also O(n) input space because it contains integers of size n 
         O(n) for the timePreferences list

    Auxiliary Space complexity:
        - The function creates two lists:
            - remaining_students: In best case the size of this list would be less than n but in worst case it can store up to n integers (student indices).
            - decreased_Preferences_time: stores up to n time preferences of the students of the 20 timeslots.
    The total space complexity is O(n)
    
          """
    remaining_students = []
    decreased_Preferences_time = []
    for i in range(n): # iterates over students 
        if preprocessed[i] == 0: # if student is already not assigned then append otherwise dont 
            remaining_students.append(i) 
            decreased_Preferences_time.append(timePreferences[i]) # append their timepreferenes too
    return remaining_students, decreased_Preferences_time

def build_reduced_proposedClasses(m, proposedClasses, pre_processed_students):
    """
    Function description:
        Constructs a modified version of the proposedClasses list that contains the reduced 
        demand in each class after pre-processing. The minimum and maximum values are updated 
        based on how many students were already assigned to each class.

    Approach description:
        For each class count how many students have already been pre-assigned.
        Subtract this count from both the original minimum and maximum constraints.
        Ensure the resulting values are not negative.
        These updated constraints are used when building the flow network.

    :Input:
        m : int
            Total number of classes.
        proposedClasses is the list of m lists where each proposedClasses[i] represnets a list with [timeslot, minimum students to be allocated, maximum students to be allocated]
        pre_processed_students : List of list
            List of (student, time_slot, class) assignments made during preprocessing.

    :returns:
        List[List[int]]
            Updated proposedClasses list with decreased min and max bounds.

    :Time complexity:
        O(m × n) where n is the number of preprocessed students.

    :Time complexity analysis:
        For each of the m classes, we iterate through all preprocessed students O(n)
        to count how many are assigned to that class. So total operations = O(m × n).
        we will remove m from the time complexity because we are not allowed to make any assumptions 
        about whether m is greater then or smaller then n according to the assignment specification 
        hence will collapse O(m + n) to O(n^2)

    :Space complexity:
        O(m)

    :Space complexity analysis:

         Input Space complexity:
        -proposedClasses has a size of O(m)
        -pre_processed_students has a size of n where n is number of preprocessed students

    Auxiliary space complexity: 

      we create a list of size m to store the classes hence the space complexity of that would be O(m)
   """
    decreased_classes = []
    for p in range(m): # Iterate over each class
        track = 0  # Counter to track how many students were pre-assigned to class p
        for s, t, c in pre_processed_students:
            if c == p:
                track = track + 1 # Count how many pre-processed students are assigned to class p
                
        # Subtract the number of pre-assigned students from the original minimum and maximum capacities        
        min_c = proposedClasses[p][1] - track
        max_c = proposedClasses[p][2] - track
        if min_c < 0:
            min_c = 0 # making sure no negatives 
        if max_c < 0:
            max_c = 0
        decreased_classes.append([proposedClasses[p][0], min_c, max_c])
    return decreased_classes

def crowdedCampus(n, m, timePreferences, proposedClasses, minimumSatisfaction):
    
    """
   Function description:
    This function aims to solve the student-to-class allocation problem under several constraints:
    - Each student must be assigned to exactly one class.
    - Each class has a minimum and maximum capacity.
    - At least `minimumSatisfaction` students must be assigned to classes whose time slots 
      fall within their top-5 preferences.

    The allocation is performed in two phases
    Firstly we perform a preprocessing phase for partial assignment and then we create a Flow network 
    through which we create a Residual network and then we run the ford Fulkerson algorithim to complete the allocation.
    Hence we can say that is uses flow network representation and the Ford-Fulkerson algorithm to find a feasible allocation that meets all constraints.

   Approach description:
    Firstly we begin with Preprocessing with pre_assign_students() function:
        _ This function attempts to satisfy the minimum capacity requirement for each class by
          assigning students whose top-5 time preferences match class time slots.
    Secondly we get unassigned students using collect_remaining_students() function:
        - Filters out students who were not pre-assigned
    Thirdly we'll check if all students were pre-assigned and satisfaction is already met, return early.

    Fourthly we'll reduce class capacities using build_reduced_proposedClasses() function:
        - This function updates each class’s capacity by subtracting the number of pre-assigned students.
        -  This helps me in ensuring that there are correct upper and lower bounds present for the remaining allocation.

    Fifthly we'll create a flow network calling Flow_Graph class:
          -In this class  I create n number of vertices for students where n represents the number of students and I create m number of vertices for class as m is the number of class then I also create 20 timeslots vertices (which are bounded). I also create a source node with negative demand to connect that with students in order to ensure that the flow is given to the students. The flow graph is initialized with reduced class capacities and remaining students because of the preprocessing now the flow graph will be created only for the remaining students in which they would be connected to 20 timeslots which are always bounded by constant which are further connected with classes with reduced capacities due to pre assigned students. These students would also be connected with source to ensure they receive the flow.

   Sixthly we'll add pre-assigned edges for the pre assigned students with fixed flow into the flow graph because they are already assigned:
        - we'll also update the demand values to reflect these fixed flows.

   Seventhly Eliminate lower bounds using eliminate_lower_bound_and_adjust_demand():
        - This step helps transform the graph into a standard flow network by adjusting demand and capacity
          so Ford-Fulkerson can operate without explicitly handling lower bounds. After eliminating the lower bounds and capacity we add edges without the lower bound between all the vertices and create super source and super sink to maintain  standard flow network. 

   Eightly we'll run ford_fulkerson():
        - In ford_fulkerson we'll further create a residual graph that we will use to augment the path.
        - we'll use a BFs  to find the augmented paths from the super source to the sink and its determines the minimum capacity of each of these paths 
        - Augments flow along the path and updates both forward and reverse edges
        - The minimum capacity of each path is the minimum flow that you can add to the flow network to augment it. 
        - Then we will also update the flow of the edges in both the residual graph and flow network by referencing 
        - It stops iterating when BFS cannot produce any more augmenting paths Stops
        - If all node demands are satisfied and total flow equals the number of unassigned students returns the updated Flow_Graph
        - Otherwise it will returns None meaning no feasible allocations are there for the remaining students 
    Ninethly: 
          we'll create a final list to show which class has a particular studnet been allocated to we'll do this by combining 
          students who were already pressigned and the students assigned through ford fulkerson.
          we'll create a list students_allocated of length n where each index will store the class number assigned to that student and then we'll first fill    in the list with pre-assigned students then we'll fill in the list for students who were assigned through the flow algorithm
        tenthly Now we'll validate the allocation:
        First by checking if every student is assigned if there is a none present if none is presnet then the allocation list is not valid so we return none  
        Then we'll how many students got their top 5 preferred time slot if the number of these students is less than minimum satisfication then we also return None
        Then it checks whether the allocation fullfills the minimum and maximum capacity if it doesn't then return None
        If either condition fails, return `None` otherwise if they fulfill all the conditions then the allocations are valid

   
   Input:
    n (int): Number of students to assign to classes.
    m (int): Number of proposed classes.
    timePreferences is the list of n lists where timePreferences[i] contains a permutation of the elements of set {0,1,...,19} to indicate the time slot preferences of student i
            proposedClasses is the list of m lists where each proposedClasses[i] represnets a list with [timeslot, minimum students to be allocated, maximum students to be allocated]
    minimumSatisfaction (int): Minimum number of students that must receive a class 
        within their top-5 preferred time slots.

  Return:
   List[int]: A list of size n where each element tells which class the student is assigned to.
   
   Time complexity:
        worst case O(n^2)
        here n is the number of students
   Time complexity Analysis:
        here n is the number of students and m is the number of classes  
       -Preprocessing students by calling the function pre_assign_students() take O(n * m) time 
       - then gathering unassigned students through the function collect_remaining_students takes 
         O(n) time
       - Reducing class capacities through build_reduced_proposedClasses takes O(n × m) time
       - Building the FlowGraph with vertices and edges by calling Flow_Graph class takes O(n + m) time 
       - Eliminating lower bounds and demands by calling function eliminate_lower_bound_and_adjust_demand also takes O(n + m) time 
       - Ford Fulkerson algorithim takes O(n × (n + m))
       - Time complexity to check if minimum statisfaction is met if all students are preprocessed is 
         O(n) because we iterate through every student to check 
       -  the loop to to add edges with fixed flow for all pre assigned students take O(n) time because it iterates 
          over the assigned students 
       - creating the students_allocated list takes O(n) because we iterate over pre_processed_students
          and Fill in the pre-assigned students in the students_allocated list and then we Fill in remaining students from flow network
          thus leading to O(n)
       - loop to check how many students were allocated to their top 5 preferences take O(n) time because we iterate over each student 
       - loop to check if the class meets the contraints is O(n + m) because first loop iterates over students 
        to count how many were assigned to that particular class and then seond loop checks if it meets the constrint
        these are two seperate loops not nested hence the time complexity is O(n + m)
        
        The total time complexity is 
          O(n x m) + O(n) +  O(n × m) + O(n + m) + O(n + m) + O(n × (n + m)) + O(n) + O(n) + O(n) + O(n) + O(n) + O(n + m)
       which will be simplified to O(n × (n + m))  + O(n x m) +  O(n × m) which will be further simplified to O(n^2)
       we will remove m from the time complexity because we are not allowed to make any assumptions 
       about whether m is greater then or smaller then n according to the assignment specification 
       hence will collapse  O(n × (n + m))  + O(n x m) +  O(n × m) to O(n^2)

  Space complexity:
       worst case: O(n)
       here n is the number of students
  
    Space complexity analysis:
        here n is the number of students and m is the number of classes
        
       Input Space
       O(n) for the timePreferences list
       O(m) for the proposedClasses list
       
       Total Input space complexity is O(n + m)
       
     Auxiliary Space complexity:
          - pre_assign_students function has total space complexity of  O(m + n)
          - collect_remaining_students function has total space complexity of O(n)
          - build_reduced_proposed_classes function has total space complexity of O(m)
          - eliminate_lower_bound_and_adjust_demand function has  O(m + n)
          - the total space complexity to build a flow graph by calling Flow_graph class is O(n + m)
          - Ford fulkerson algorithim takes O(n + m) auxiliary space 
          - we create allocate array with size O(n) 
          - students_allocated array also takes up size O(n) 
          - class_count array takes O(m) size 
          
          hence the total space complexity is O(n + m) + O(n + m) + O(n + m) + O(n + m) + O(n + m) +O(n) + O(m) + O(n) + O(n) + O(m)
          which we can simplify to O(n + m) which would further simplify to O(n)
          as we will remove m from the time complexity because we are not allowed to make any assumptions 
          about whether m is greater then or smaller then n according to the assignment specification 
          hence will collapse O(n + m) to O(n) 

"""
    # firstly we'll preprocess to pre-assign students to help to their top 5 preferences 
    pre_processed_students, preprocessed, num_satisfied = pre_assign_students(n, m, timePreferences, proposedClasses)
    
    # get students who did not get assigned during pre-processing
    remaining_students, decreased_Preferences_time = collect_remaining_students(n, preprocessed, timePreferences)
    number_remain_students = len(remaining_students)

    # If all students have been preprocessed
    if number_remain_students == 0:
        # Check if minimum satisfaction is achieved
        if num_satisfied >= minimumSatisfaction:
            allocate = [None] * n
            for s, t, c in pre_processed_students:
                allocate[s] = c
            return allocate
        else:
            return None # Not enough students got their top 5 timeslots so terminate
        
    # Now Update class capacities by subtracting students already assigned earlier
    decreased_classes = build_reduced_proposedClasses(m, proposedClasses, pre_processed_students)
   
    # check how many more satisfied students are needed
    satisfaction_achieved = minimumSatisfaction - num_satisfied
    if satisfaction_achieved < 0:
        satisfaction_achieved = 0
   
    # now create the flow network using remaining students
    create_flow_network = Flow_Graph(number_remain_students, m, decreased_Preferences_time, decreased_classes)

    # Add edges for the students already pre-assigned in preprocessing
    for s, t, c in pre_processed_students:
        student_ver = create_flow_network.get_vertex(s + 1)
        timeslot_ver = create_flow_network.get_vertex(create_flow_network.n + t + 1)
        class_ver = create_flow_network.get_vertex(create_flow_network.n + 20 + c + 1)

        create_flow_network.source_node.edges.append(Flow_Edge(student_ver, 1, 1))
        student_ver.edges.append(Flow_Edge(timeslot_ver, 1, 1)) # add edges with fixed flow
        timeslot_ver.edges.append(Flow_Edge(class_ver, 1, 1))
        
        # Update demand to show pre-assignment
        create_flow_network.demand[0] = create_flow_network.demand[0] - 1
        create_flow_network.demand[class_ver.id] = create_flow_network.demand[class_ver.id] + 1

    # Eliminate lower bounds and adjust demands before running the Ford Fulkerson algorithim
    create_flow_network.eliminate_lower_bound_and_adjust_demand()
    create_ford_flow_network = ford_fulkerson(create_flow_network) # Run Ford Fulkerson to allocate remaining students
    if create_ford_flow_network is None:
        return None # No feasible solution
  
    students_allocated = [None] * n
    for s, t, c in pre_processed_students: # Fill in the pre-assigned students in the allocated list to show allocation
        students_allocated[s] = c

    for s in range(number_remain_students): # Fill in remaining students from flow network
        student = create_ford_flow_network.student_nodes[s][0]
        for stu_e in student.edges:
            if stu_e.flow == 1:
                time_s = stu_e.v
                for time_ed in time_s.edges:
                    if time_ed.flow == 1:
                        c = time_ed.v.id - (number_remain_students + 20 + 1)
                        stu_id = remaining_students[s]
                        students_allocated[stu_id] = c
                        break
                break

    

    # counting the amount of students who get their top 5 preferences of the timeslot
    satisfied = 0
    for stu in range(n):
        assigned_class = students_allocated[stu]
        if assigned_class is not None:
            class_time = proposedClasses[assigned_class][0]
            i = 0
            while i < 5 and i < len(timePreferences[stu]):
                if timePreferences[stu][i] == class_time:
                    satisfied = satisfied + 1
                    break
                i = i + 1
    if satisfied < minimumSatisfaction: # If not enough students are satisfied return none
        return None
    
    if None in students_allocated or len(students_allocated) == 0: # validating if  every student is assigned
        return None
    
    class_counts = [0] * m
    for assigned_class in students_allocated: # loop to check how many students are assigned to this class 
        if assigned_class is not None:  
            class_counts[assigned_class] = class_counts[assigned_class] + 1

    for class_id in range(m): # loop to check if the class meets the constraints
        min_cap, max_cap = proposedClasses[class_id][1], proposedClasses[class_id][2]
        if class_counts[class_id] < min_cap or class_counts[class_id] > max_cap:
            return None

    return students_allocated

# QUESTION 2




class Node:
    """
    class Description:
        This class represents a node in a trie structure where each character of a word is stored.
    """
    def __init__(self,payload=None,size=27):
        """
    Function description:
    Initializes a trie node which represents a character of a word for the Bad_AI class. Each node has a fixed number of attributes like 
    a list of fixed size of child nodes, payload for storing character and path length, and a field for storing a entire word ( this also acts as a payload)

    
    :Input:
    payload:tuple of (char, path_len): representing data stored at each node in for trie. which helps during traversal
    size:(int): this is an integer on which the size of the link list is built. 
                    The size is 27 to represent alphabets from a to z and one special character for the terminal
    postcondition:
    A new Node instance is created with with its attributes

    :Time complexity:
    O(1)

    :Time complexity analysis:

    Time complexity is O(1) because the constructor performs constant operation of assigning values

    :Space complexity:
    O(1)

    :Space complexity analysis:
    Input space complexity:
       payload stores a tuple which takes constant space O(1)
       and size is a single value integer so also O(1)
    Auxilairy Space complexity:

    Each node uses a constant has a list of fixed size and two other fields hence we can say that the space is O(1)
    """
        self.link = [None] * size
        self.payload = None   # (char, path_len) tuple for non-terminal nodes
        self.word_here = None # word at terminal node

class Trie:
    """
    Class description:
        
    This class creates a trie structure which is used for storing words. so that later Bad_AI class can make use of this class to look for words that differ 
    by 1 substitution. This class manages a prefix tree where each path from the root node to a terminal node 
    represents a complete word.
    I chose Trie as a data structure for this question because it enables you to search character by character at a time which makes it perfect for this question because we can detect substitutions of characters at positions
    """
    def __init__(self):
        """
        Function description:
        Initializes the Trie by creating its root node. This root node serves as the starting
        point for inserting and searching all words in the trie.
    
        :Input:
        None
    
        :postcondition:
        A Trie object with one root node is initialized.
    
        :Time complexity:
        O(1)
    
        :Time complexity analysis:
        The constructor performs a single Node creation which is a constant operation of O(1)
    
        :Space complexity:
        O(1)
    
        :Space complexity analysis:
        Only one Node object is created hence the space required is also constant of O(1) as it requires fixed amount of space.
      """
        self.start = Node()

    def insert(self, word):
        """
    Function description:
    Inserts a word into the trie by creating nodes for each character. It also records
    payload at each node and marks the end of the word using a special terminal node.

    Approach description (if main function):
    Starting from the root node each character of the word is mapped to an index has an index between 1 to 26,
    and the corresponding child node is created if it doesn't exist. Each node's payload stores
    the character and the number of characters remaining to the end of the word.
    After all characters are inserted a terminal node is inserted at the index 0 this is done to ensure the in order traversal. The terminal node also
    contains a special kind of payload which is the entire word from the root to that particular terminal.
    :Input:(string)
    word: A lowercase string consisting of letters a–z that has to be inserted into the trie.

    :postcondition:
    The word is stored in the trie. 

    :Time complexity:
            Worst case: O(M), where M is the number of character in the longest word we added (length of the word)
        
    :Time complexity analysis:
    - The loop runs M times once for each character.
    - In each iteration we do constant-time operations like converting the character to an index, creating a node for the character and setting the payload at that particular node like char and path length which is the number of characters remaining in the word after the current character. Adding a terminal node and storing the entire word at that terminal node also takes O(1).
      Hence the total time complexity is O(M)

    :Space complexity:
       worse case O(M)

    :Space complexity analysis:
     
     Input Space complexity:
      The input word is a string of M lowercase letters (a–z) hence it will require O(M) space.
     Auxiliary Space Complexity:
    
    - In the worst case there will be no shared prefixes hence M new nodes will be created for the characters
      plus one terminal node at the end.
    - Each node has a fixed array of 27 links and a few attributes so space per node is constant.
     
      Hence the total space complexity is O(M)
     """
        node = self.start
        L = len(word)
        for i in range(L):
            index = ord(word[i]) - 97 + 1
            if node.link[index] is None: 
                node.link[index] = Node() # Create new node if it doesn't exist
            node = node.link[index]
            
            node.payload = (word[i], L - i - 1) # Store (char, remaining chars to end)
        if node.link[0] is None:
            node.link[0] = Node() # Add terminal node at link[0] to mark end of word
        node.link[0].word_here = word  # Store full word at terminal
        
        
class Bad_AI:
    """
    class description:
        This class represents the output of a faulty AI which gives incorrect words that are
        exactly one subsitution away from a the given suspicious word. This calls helps identify words
        that differ by only one character from the given word and it only allows subsitution no insertion
        or deletion is allowed
    """
    def __init__(self, wordlist):
        """
    Function description:
        This function creates a trie data structure and then inserts every word from list_words
        into the trie character-by-character .

    Approach description:
        This function is the constructor of the Bad_AI class. It initializes an empty Trie
        and then iterates over each word in list_words. Each word is inserted into the trie
        where every character of the word is added as a node in the trie. A
        terminal node marks the end of a word. 

    Input:
        list_words (list): A list of lowercase words that need to be stored in the trie data structure.

    PostCondition:
        All the words are to be inserted into the trie.

    Time complexity:
        Worst case complexity: O(C)
        where C is the total number of characters across all words in list_words.

    Time complexity analysis:
        The function iterates over N words in list_words and each word has a length of M hence when it performs
        an insert operation into the trie it takes O(M) time since each character is processed one at a time. Hence for N words this insert opertion
        will takes O(M * N) time. We can also say that all the characters of all the words add up to hence we have C characters in total to be added to the     Trie. Hence in worst case scenerio we will create C nodes in the Trie thus the time complexity will be O(C) as So O(C) <= O(M*N)
    Space Complexity:
           worst case O(C) 
          
     where C is the total number of characters across all words in list_words.
     
    Space complexity Analysis:
           
  
    Input space: O(C) where C is the number of characters across all input words.

    The auxiliary space: is due to constructing the trie charachter by character and terminal node is stored as a
    separate node and in the worst case where no
    words share a prefix O(C) nodes are created. Hence total auxiliary space complexity
    is O(C).
    
     """
        self.trie = Trie()
        for w in wordlist:
            self.trie.insert(w)  # Insert each word into trie

    def check_word(self, sus_word):
        """
        Function description:
        This function is responsible for finding and returning all words in the trie that differ from
        the given sus_word by exactly one character substitution. It does so by searching the trie
        in an in order way using a stack-based array to keep track of the current traversal state.
    
        During search the function maintains three pieces of state:
         The current trie node being visited.
         The current position of index in the suspect word.
         The number of mismatched characters so far.
    
        The function allows only one difference of characters in words between the two words (substitution) and prunes branches that have more than 1 difference
        If a path reaches a terminal node with exactly one substitution (1 difference) and the suspect word has been has reached its last character
        Then this word is considered valid and then added to the result list.
        
        Approach Description:
    The function performs an in-order search of the trie using a stack-based array called discovered. As stack helps us explore all the full word in depth like exploring child node of child node and the proceeding to the next branch this way we can complete our entire word in the trie by going deep.
    Each element in the stack is a tuple consisting of:
      - The current trie node being visited.
      - The current character position in the suspect word.
      - The number of mismatched characters.

    The search starts by appending the initial state (root node, index 0, 0 mismatches) to the stack.
    While the stack is not empty the function repeatedly pops the most recent state (LIFO order) and
    explores the child nodes of the current trie node from index 1 to 26, which are alphabets a to z ( this ensures alphabetical order). Index 0 is reserved for the terminal node which holds complete words.

    For each valid child node of the current node:
    - The function retrieves the stored character and path length from its payload.
    - If the character from the payload of the child node matches of the character at the sus_word at the current position the search continues and the this child node is appended into the stack with the position of the current character of the sus-word and number of mismatches.
    - If the character does not match but no previous mismatch has occurred
      the wrong count is incremented to 1 and the child node is also append.
    - If a mismatch has already occurred wrong == 1 then the function only continues if the
      remaining characters in sus_word can still be matched within the
      remaining trie path length path_len. If not then we move to the next branch.

    Once the length of the sus_word becomes 0 the function checks whether
    the current node has a terminal child at index 0. If it does and exactly one there is only one mismatch between the charactes of the words then 
    the stored word in word_here is appended to the result list otherwise if no word has 1 subsitution form the sus_word then it doesn't append anything and the result list remains empty 
            
    input 
    sus word (string) : word we need to compare the words in trie against
    
    Time complexity:
        worstcase O(J ∗ N)+O(X) 
        where J is the length of the sus word and N is the number of words in the Trie and X is the number of characters in the result list 
    
    Time complexity analysis:
       In the worst case for each of the J characters in sus_word, the function may search upto
       26 children at each level of the trie. Since the trie could store up to N different words,
       Thus this would make the complexity the worst-case O(J × N).
       o(X) will be the complexity of creating the result list which contains X number of characters
       
       Total time complexity is   O(J × N) + O(X)
       
     space complexity:worse case O(X)
     where J is the length of the sus word and X is the number of characters in the result list 
     
     Input Space complexity is O(1) because sus_word is a single string of O(1)
     
     
     Auxiliary space analysis:
     
     The discovered stack is used to manage in-order search. Each element in the stack is a fixed-size
       tuple consisting of the current node, the index of the current character of the sus_word, and the mismatch count. Since the stack
       depth is bounded by the length of sus_word (J) its size in the worst case scenario will also be O(J) which remains small and
       does not grow with the number of matches. The result list stores words that match that have a substitution of 1.
       If there are many matches this list can store up to X characters in total leading to O(X) auxiliary space. Since O(X) dominates O(J) in the worst case the overall auxiliary space complexity is O(X).
    
    Total Space complexity is O(X)
    """
        result = []
        discovered = [(self.trie.start, 0, 0)]
     
        while len(discovered) > 0:
            node, j, wrong = discovered.pop()
            chars_left = len(sus_word) - j
     
            if chars_left == 0: #reached end of sus_word
                last_node = node.link[0] # checking terminal 
                if last_node and last_node.word_here and wrong == 1:
                    result.append(last_node.word_here) # Valid match
                continue
     
            for c in range(1, 27):
                kid = node.link[c]
                if not kid or not kid.payload: # invalid child
                    continue
     
                node_char, node_len = kid.payload
     
                # if there are more then 1 mismatch and the path is too short move to the next branch
                if wrong == 1 and node_len < chars_left - 1:
                    continue
     
                # Compare character from payload 
                if sus_word[j] == node_char: # If no mismatch
                    discovered.append((kid, j + 1, wrong))
                elif wrong == 0:  # allow one mismatch
                    discovered.append((kid, j + 1, 1))
        return result

    
 













